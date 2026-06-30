"""
api.py
------
FastAPI REST API for the Drug Target Explorer.

Run:
    uvicorn src.api:app --reload --port 8000

Endpoints
    GET /health
    GET /targets
    GET /drugs/{target}
    GET /drugs/{target}/top
"""

import sqlite3
from pathlib import Path
from typing import Optional

import pandas as pd
from fastapi import FastAPI, HTTPException, Query
from pydantic import BaseModel

DB_PATH = Path(__file__).resolve().parent.parent / "drug_targets.db"

app = FastAPI(
    title="Drug Target Explorer",
    description="Return drugs for a target, ranked by specificity and adverse-effect burden.",
    version="1.0.0",
)


def get_conn():
    if not DB_PATH.exists():
        raise HTTPException(status_code=503, detail="Database not found. Run `python src/etl.py` first.")
    return sqlite3.connect(DB_PATH)


class DrugResult(BaseModel):
    rank: int
    drug_id: str
    drug_name: str
    n_targets: int
    median_llr: Optional[float]
    composite_score: float
    is_approved: Optional[int] = None
    max_clinical_phase: Optional[float] = None
    black_box_warning: Optional[int] = None
    has_been_withdrawn: Optional[int] = None


class TargetDrugsResponse(BaseModel):
    target_query: str
    resolved_ensembl_id: Optional[str]
    resolved_symbol: Optional[str]
    n_drugs_found: int
    drugs: list[DrugResult]


def resolve_target(raw: str, con) -> tuple[str, str]:
    df = pd.read_sql(
        """
        SELECT DISTINCT ensembl_id, gene_symbol
        FROM drug_targets
        WHERE ensembl_id = ?
           OR UPPER(gene_symbol) = UPPER(?)
        """,
        con,
        params=(raw, raw),
    )

    if df.empty:
        raise HTTPException(
            status_code=404,
            detail=f"Target '{raw}' not found. Use /targets to inspect available targets.",
        )

    return df.iloc[0]["ensembl_id"], df.iloc[0]["gene_symbol"]


def score_drugs_for_target(ensembl_id: str, con) -> pd.DataFrame:
    target_drugs = pd.read_sql(
        """
        SELECT DISTINCT drug_id
        FROM drug_targets
        WHERE ensembl_id = ?
        """,
        con,
        params=(ensembl_id,),
    )

    if target_drugs.empty:
        return pd.DataFrame()

    ids = target_drugs["drug_id"].tolist()
    placeholders = ",".join("?" * len(ids))

    names_and_meta = pd.read_sql(
        f"""
        SELECT
            drug_id,
            drug_name,
            is_approved,
            max_clinical_phase,
            black_box_warning,
            has_been_withdrawn
        FROM drugs
        WHERE drug_id IN ({placeholders})
        """,
        con,
        params=ids,
    )

    specificity = pd.read_sql(
        f"""
        SELECT
            drug_id,
            COUNT(DISTINCT ensembl_id) AS n_targets
        FROM drug_targets
        WHERE drug_id IN ({placeholders})
        GROUP BY drug_id
        """,
        con,
        params=ids,
    )

    ae = pd.read_sql(
        f"""
        SELECT drug_id, llr
        FROM adverse_effects
        WHERE drug_id IN ({placeholders})
        """,
        con,
        params=ids,
    )

    if ae.empty:
        ae_agg = pd.DataFrame(columns=["drug_id", "median_llr"])
    else:
        ae["llr"] = pd.to_numeric(ae["llr"], errors="coerce")
        ae_agg = (
            ae.groupby("drug_id", as_index=False)["llr"]
            .median()
            .rename(columns={"llr": "median_llr"})
        )

    df = (
        target_drugs
        .merge(names_and_meta, on="drug_id", how="left")
        .merge(specificity, on="drug_id", how="left")
        .merge(ae_agg, on="drug_id", how="left")
    )

    df["drug_name"] = df["drug_name"].fillna(df["drug_id"])
    df["n_targets"] = df["n_targets"].fillna(9999).astype(int)
    df["median_llr"] = pd.to_numeric(df["median_llr"], errors="coerce")

    df["rank_specificity"] = df["n_targets"].rank(method="average", ascending=True)
    df["rank_safety"] = df["median_llr"].rank(method="average", ascending=True, na_option="bottom")

    approval_bonus = df["is_approved"].fillna(0).astype(float) * -0.25
    phase_bonus = df["max_clinical_phase"].fillna(0).astype(float) * -0.05
    warning_penalty = df["black_box_warning"].fillna(0).astype(float) * 0.25
    withdrawn_penalty = df["has_been_withdrawn"].fillna(0).astype(float) * 1.0

    df["composite_score"] = (
        df["rank_specificity"]
        + df["rank_safety"]
        + warning_penalty
        + withdrawn_penalty
        + approval_bonus
        + phase_bonus
    )

    df = df.sort_values(
        by=["composite_score", "n_targets", "median_llr"],
        ascending=[True, True, True],
    ).reset_index(drop=True)

    df["rank"] = df.index + 1
    return df


@app.get("/health")
def health():
    return {"status": "ok", "db_path": str(DB_PATH)}


@app.get("/targets")
def list_targets(limit: int = Query(default=100, ge=1, le=5000)):
    con = get_conn()
    df = pd.read_sql(
        """
        SELECT
            ensembl_id,
            gene_symbol,
            COUNT(DISTINCT drug_id) AS n_drugs
        FROM drug_targets
        GROUP BY ensembl_id, gene_symbol
        ORDER BY n_drugs DESC, gene_symbol ASC
        LIMIT ?
        """,
        con,
        params=(limit,),
    )
    con.close()
    return df.to_dict(orient="records")


@app.get("/drugs/{target}", response_model=TargetDrugsResponse)
def get_drugs_for_target(
    target: str,
    top_n: int = Query(default=20, ge=1, le=500),
):
    con = get_conn()
    ensembl_id, gene_symbol = resolve_target(target, con)
    df = score_drugs_for_target(ensembl_id, con)
    con.close()

    if df.empty:
        return TargetDrugsResponse(
            target_query=target,
            resolved_ensembl_id=ensembl_id,
            resolved_symbol=gene_symbol,
            n_drugs_found=0,
            drugs=[],
        )

    df = df.head(top_n)

    drugs = []
    for row in df.itertuples():
        drugs.append(
            DrugResult(
                rank=int(row.rank),
                drug_id=str(row.drug_id),
                drug_name=str(row.drug_name),
                n_targets=int(row.n_targets),
                median_llr=None if pd.isna(row.median_llr) else round(float(row.median_llr), 4),
                composite_score=round(float(row.composite_score), 4),
                is_approved=None if pd.isna(row.is_approved) else int(row.is_approved),
                max_clinical_phase=None if pd.isna(row.max_clinical_phase) else float(row.max_clinical_phase),
                black_box_warning=None if pd.isna(row.black_box_warning) else int(row.black_box_warning),
                has_been_withdrawn=None if pd.isna(row.has_been_withdrawn) else int(row.has_been_withdrawn),
            )
        )

    return TargetDrugsResponse(
        target_query=target,
        resolved_ensembl_id=ensembl_id,
        resolved_symbol=gene_symbol,
        n_drugs_found=len(drugs),
        drugs=drugs,
    )


@app.get("/drugs/{target}/top", response_model=DrugResult)
def get_top_drug_for_target(target: str):
    response = get_drugs_for_target(target=target, top_n=1)
    if not response.drugs:
        raise HTTPException(status_code=404, detail=f"No drugs found for target '{target}'.")
    return response.drugs[0]