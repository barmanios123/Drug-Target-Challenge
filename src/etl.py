"""
etl.py
------
Reads data/molecules.tsv and data/adverseEffects.tsv,
resolves Ensembl IDs to HGNC symbols via biothings_client,
and writes normalized tables into drug_targets.db (SQLite).

This version is tuned for the Pfizer coding challenge dataset:
- uses molecules.id as the primary drug identifier
- parses linkedtargets safely from bracketed list-like strings
- normalises adverseEffects drug IDs to CHEMBL-prefixed format
- keeps drug metadata useful for downstream API queries

Usage
    python src/etl.py
"""

import logging
import sqlite3
from pathlib import Path

import pandas as pd
import biothings_client as bt

logging.basicConfig(level=logging.INFO, format="%(levelname)s  %(message)s")
log = logging.getLogger(__name__)

ROOT    = Path(__file__).resolve().parent.parent
DATA    = ROOT / "data"
DB_PATH = ROOT / "drug_targets.db"


def read_tsv(name: str) -> pd.DataFrame:
    for suffix in (name, name + ".gz"):
        p = DATA / suffix
        if p.exists():
            log.info("Reading %s", p.name)
            return pd.read_csv(
                p,
                sep="\t",
                low_memory=False,
                compression="gzip" if suffix.endswith(".gz") else None,
            )
    raise FileNotFoundError(f"'{name}' not found in {DATA}")


def normalise_cols(df: pd.DataFrame) -> pd.DataFrame:
    df.columns = (
        df.columns.astype(str)
        .str.strip()
        .str.lower()
        .str.replace(r"[\s\-/]+", "_", regex=True)
    )
    return df


def resolve_ensembl(ids: list[str]) -> dict[str, str]:
    clean = sorted({str(i).split(".")[0] for i in ids if str(i).startswith("ENSG")})
    log.info("Resolving %d Ensembl IDs via biothings...", len(clean))

    mg = bt.get_client("gene")
    results = mg.querymany(
        clean,
        scopes="ensembl.gene",
        fields="symbol",
        species="human",
    )

    mapping = {}
    for r in results:
        query  = r.get("query")
        symbol = r.get("symbol")
        if query and symbol:
            mapping[query] = symbol

    n_missing = len(clean) - len(mapping)
    if n_missing:
        log.warning("%d Ensembl IDs were not resolved to symbols", n_missing)

    return mapping


def parse_target_list(raw) -> list[str]:
    if pd.isna(raw):
        return []

    cleaned = (
        str(raw)
        .replace("[", "")
        .replace("]", "")
        .replace("'", "")
        .replace('"', "")
    )

    targets = []
    for token in cleaned.replace(",", "|").split("|"):
        token = token.strip()
        if token.startswith("ENSG"):
            targets.append(token.split(".")[0])

    return list(dict.fromkeys(targets))


def parse_molecules(df: pd.DataFrame):
    df = normalise_cols(df)
    log.info("molecules columns: %s", df.columns.tolist())

    id_col     = "id" if "id" in df.columns else df.columns[0]
    name_col   = "name" if "name" in df.columns else None
    target_col = "linkedtargets" if "linkedtargets" in df.columns else None

    if target_col is None:
        for c in df.columns:
            sample = df[c].dropna().astype(str).head(50)
            if sample.str.contains(r"ENSG\d{11}", regex=True).any():
                target_col = c
                break

    if target_col is None:
        raise ValueError("Could not find target column in molecules.tsv")

    log.info("id=%s | name=%s | targets=%s", id_col, name_col, target_col)

    drug_cols    = [id_col]
    optional_cols = [
        "name",
        "drugtype",
        "maximumclinicaltrialphase",
        "isapproved",
        "hasbeenwithdrawn",
        "blackboxwarning",
        "yearoffirstapproval",
    ]
    for c in optional_cols:
        if c in df.columns and c not in drug_cols:
            drug_cols.append(c)

    drugs = df[drug_cols].drop_duplicates().rename(
        columns={
            id_col: "drug_id",
            "name": "drug_name",
            "drugtype": "drug_type",
            "maximumclinicaltrialphase": "max_clinical_phase",
            "isapproved": "is_approved",
            "hasbeenwithdrawn": "has_been_withdrawn",
            "blackboxwarning": "black_box_warning",
            "yearoffirstapproval": "year_of_first_approval",
        }
    )

    if "drug_name" not in drugs.columns:
        drugs["drug_name"] = drugs["drug_id"]

    rows = []
    for _, row in df.iterrows():
        drug_id = str(row[id_col])
        targets = parse_target_list(row[target_col])
        for t in targets:
            rows.append({"drug_id": drug_id, "ensembl_id": t})

    dt = pd.DataFrame(rows).drop_duplicates()

    if dt.empty:
        raise ValueError("No drug-target links were parsed from molecules.tsv")

    return drugs, dt


def parse_adverse_effects(df: pd.DataFrame) -> pd.DataFrame:
    df = normalise_cols(df)
    log.info("adverseEffects columns: %s", df.columns.tolist())

    id_col    = "chembl_id" if "chembl_id" in df.columns else df.columns[0]
    event_col = "event" if "event" in df.columns else None
    llr_col   = "llr"   if "llr"   in df.columns else None

    if event_col is None:
        event_col = next(
            (c for c in df.columns if any(k in c for k in ("event", "effect", "adr", "adverse"))),
            None,
        )

    if llr_col is None:
        llr_col = next(
            (c for c in df.columns if "llr" in c or "score" in c or "ratio" in c),
            None,
        )

    log.info("id=%s | event=%s | llr=%s", id_col, event_col, llr_col)

    out = df.rename(columns={id_col: "drug_id"})

    if event_col:
        out = out.rename(columns={event_col: "event"})
    else:
        out["event"] = "unknown"

    if llr_col:
        out = out.rename(columns={llr_col: "llr"})
    else:
        out["llr"] = float("nan")

    # --- KEY FIX: normalise AE drug IDs to CHEMBL-prefixed format ---
    out["drug_id"] = out["drug_id"].astype(str).str.strip()
    out["drug_id"] = out["drug_id"].apply(
        lambda x: x if x.upper().startswith("CHEMBL") else f"CHEMBL{x}"
    )

    out["event"] = out["event"].astype(str)
    out["llr"]   = pd.to_numeric(out["llr"], errors="coerce")

    return out[["drug_id", "event", "llr"]].dropna(subset=["drug_id"])


def write_db(drugs: pd.DataFrame, dt: pd.DataFrame, ae: pd.DataFrame):
    log.info("Writing %s ...", DB_PATH)
    con = sqlite3.connect(DB_PATH)

    drugs.to_sql("drugs",           con, if_exists="replace", index=False)
    dt.to_sql   ("drug_targets",    con, if_exists="replace", index=False)
    ae.to_sql   ("adverse_effects", con, if_exists="replace", index=False)

    cur = con.cursor()
    cur.execute("CREATE INDEX IF NOT EXISTS idx_drugs_id   ON drugs(drug_id)")
    cur.execute("CREATE INDEX IF NOT EXISTS idx_dt_ensembl ON drug_targets(ensembl_id)")
    cur.execute("CREATE INDEX IF NOT EXISTS idx_dt_symbol  ON drug_targets(gene_symbol)")
    cur.execute("CREATE INDEX IF NOT EXISTS idx_dt_drug    ON drug_targets(drug_id)")
    cur.execute("CREATE INDEX IF NOT EXISTS idx_ae_drug    ON adverse_effects(drug_id)")
    con.commit()
    con.close()

    log.info(
        "Done — %d drugs | %d drug-target links | %d AE records",
        len(drugs), len(dt), len(ae),
    )


def run():
    mol_raw = read_tsv("molecules.tsv")
    ae_raw  = read_tsv("adverseEffects.tsv")

    drugs, dt = parse_molecules(mol_raw)

    mapping       = resolve_ensembl(dt["ensembl_id"].unique().tolist())
    dt["gene_symbol"] = dt["ensembl_id"].map(mapping).fillna("UNKNOWN")

    ae = parse_adverse_effects(ae_raw)
    write_db(drugs, dt, ae)


if __name__ == "__main__":
    run()