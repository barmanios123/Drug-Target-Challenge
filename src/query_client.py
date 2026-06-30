"""
query_client.py
---------------
Bioinformatician-facing CLI script.
Queries the Drug Target Explorer API, prints a ranked table,
and saves a bubble chart PNG and ranking CSV to the output directory.

Usage
    python src/query_client.py AR
    python src/query_client.py ENSG00000169083 --top 10 --out results/
"""

import argparse
import sys
from pathlib import Path

import pandas as pd
import requests
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from rich.console import Console
from rich.table import Table

BASE_URL = "http://127.0.0.1:8000"
console  = Console()


def fetch_drugs(target: str, top_n: int) -> dict:
    url = f"{BASE_URL}/drugs/{target}?top_n={top_n}"
    try:
        r = requests.get(url, timeout=30)
    except requests.ConnectionError:
        console.print("[bold red]Cannot reach API.[/] "
                      "Start it with: uvicorn src.api:app --reload --port 8000")
        sys.exit(1)

    if r.status_code == 404:
        console.print(f"[bold red]Target '{target}' not found.[/] "
                      f"Try: curl {BASE_URL}/targets")
        sys.exit(1)

    r.raise_for_status()
    return r.json()


def bubble_chart(df: pd.DataFrame, symbol: str, out_dir: Path) -> Path:
    """
    X = n_targets       (specificity  - fewer = more specific)
    Y = median_llr      (safety       - lower = safer)
    Bubble size         = inversely proportional to composite_score
    Colour              = composite rank (green = best)
    """
    df = df.copy()

    llr_max  = df["median_llr"].max(skipna=True)
    fill_val = float(llr_max * 1.15) if pd.notna(llr_max) else 1.0
    df["median_llr_plot"] = df["median_llr"].astype(float).fillna(fill_val)

    df["bubble_size"] = (
        (1.0 / df["composite_score"].clip(lower=0.1)) * 12000
    ).clip(upper=3000)

    fig, ax = plt.subplots(figsize=(10, 6))

    sc = ax.scatter(
        df["n_targets"],
        df["median_llr_plot"],
        s=df["bubble_size"],
        c=df["rank"],
        cmap="RdYlGn_r",
        alpha=0.75,
        edgecolors="white",
        linewidths=0.6,
        zorder=3,
    )

    for _, row in df.iterrows():
        ax.annotate(
            row["drug_name"],
            (row["n_targets"], row["median_llr_plot"]),
            fontsize=6.5,
            ha="center",
            va="bottom",
            xytext=(0, 5),
            textcoords="offset points",
            color="#333333",
        )

    cbar = plt.colorbar(sc, ax=ax, pad=0.02)
    cbar.set_label("Composite rank  (lower = better)", fontsize=9)

    ax.set_xlabel("Number of targets hit  (fewer = more specific)", fontsize=10)
    ax.set_ylabel("Median LLR of adverse events  (lower = safer)", fontsize=10)
    ax.set_title(f"Drug landscape for target: {symbol}", fontsize=13, fontweight="bold")
    ax.xaxis.set_major_locator(ticker.MaxNLocator(integer=True))
    ax.grid(True, linestyle="--", linewidth=0.4, alpha=0.5, zorder=0)

    q25_x = df["n_targets"].quantile(0.25)
    q25_y = df["median_llr_plot"].quantile(0.25)
    ax.axvline(q25_x, color="steelblue", linestyle="--", linewidth=0.8, alpha=0.5)
    ax.axhline(q25_y, color="green",     linestyle="--", linewidth=0.8, alpha=0.5)
    ax.text(0.02, 0.03, "ideal zone (specific + safe)",
            transform=ax.transAxes, fontsize=8, color="gray")

    plt.tight_layout()
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"drugs_{symbol}.png"
    plt.savefig(out_path, dpi=150)
    plt.close()
    return out_path


def save_ranking_table(df: pd.DataFrame, symbol: str, out_dir: Path) -> Path:
    """Save a clean ranking CSV to the results folder."""
    out = df[[
        "rank", "drug_id", "drug_name", "n_targets",
        "median_llr", "composite_score", "is_approved",
        "max_clinical_phase", "black_box_warning", "has_been_withdrawn"
    ]].copy()

    out.columns = [
        "Rank", "Drug ID", "Drug Name", "# Targets",
        "Median LLR", "Composite Score", "Approved",
        "Max Clinical Phase", "Black Box Warning", "Withdrawn"
    ]

    out["Approved"]          = out["Approved"].apply(lambda x: "Yes" if x == 1 else "No")
    out["Black Box Warning"] = out["Black Box Warning"].apply(lambda x: "Yes" if x == 1 else "No")
    out["Withdrawn"]         = out["Withdrawn"].apply(lambda x: "Yes" if x == 1 else "No")
    out["Median LLR"]        = out["Median LLR"].apply(
        lambda x: f"{x:.4f}" if pd.notna(x) else "N/A"
    )

    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"rankings_{symbol}.csv"
    out.to_csv(out_path, index=False)
    return out_path


def print_table(data: dict) -> pd.DataFrame:
    symbol  = data.get("resolved_symbol") or data["target_query"]
    ensembl = data.get("resolved_ensembl_id", "")
    drugs   = data["drugs"]

    console.print(
        f"\n[bold cyan]Target:[/] {data['target_query']} "
        f"[bold]{symbol}[/] ({ensembl})\n"
        f"[bold cyan]{data['n_drugs_found']}[/] drugs found\n"
    )

    tbl = Table(show_header=True, header_style="bold magenta", show_lines=False)
    tbl.add_column("Rank",       justify="right",  style="bold")
    tbl.add_column("Drug ID",    justify="left")
    tbl.add_column("Drug Name",  justify="left")
    tbl.add_column("# Targets",  justify="right")
    tbl.add_column("Median LLR", justify="right")
    tbl.add_column("Approved",   justify="center")
    tbl.add_column("Phase",      justify="center")
    tbl.add_column("Score",      justify="right",  style="dim")

    for d in drugs:
        llr_str      = f"{d['median_llr']:.3f}" if d["median_llr"] is not None else "N/A"
        approved_str = "Y" if d.get("is_approved") == 1 else "-"
        phase_val    = d.get("max_clinical_phase")
        phase_str    = str(int(phase_val)) if phase_val is not None else "-"

        tbl.add_row(
            str(d["rank"]),
            d["drug_id"],
            d["drug_name"],
            str(d["n_targets"]),
            llr_str,
            approved_str,
            phase_str,
            f"{d['composite_score']:.2f}",
        )

    console.print(tbl)

    best = drugs[0]
    console.print(
        f"\n[bold green]Best drug:[/] [bold]{best['drug_name']}[/] ({best['drug_id']})\n"
        f"   Targets: {best['n_targets']}  |  "
        f"Median LLR: {best['median_llr'] if best['median_llr'] is not None else 'N/A'}  |  "
        f"Approved: {'Yes' if best.get('is_approved') == 1 else 'No'}  |  "
        f"Phase: {int(best['max_clinical_phase']) if best.get('max_clinical_phase') is not None else '-'}  |  "
        f"Score: {best['composite_score']}\n"
    )

    return pd.DataFrame(drugs)


def main():
    parser = argparse.ArgumentParser(
        description="Query Drug Target Explorer for a given gene target."
    )
    parser.add_argument("target",
        help="Ensembl ID (e.g. ENSG00000169083) or gene symbol (e.g. AR)")
    parser.add_argument("--top",  type=int, default=20,
        help="Number of top drugs to retrieve (default 20)")
    parser.add_argument("--out",  default="results",
        help="Output directory for chart PNG and ranking CSV (default: ./results/)")
    args = parser.parse_args()

    data   = fetch_drugs(args.target, args.top)
    df     = print_table(data)
    symbol = data.get("resolved_symbol") or args.target

    if not df.empty:
        out_dir    = Path(args.out)
        chart_path = bubble_chart(df, symbol, out_dir)
        csv_path   = save_ranking_table(df, symbol, out_dir)
        console.print(f"[dim]Chart saved:   {chart_path}[/]")
        console.print(f"[dim]Ranking saved: {csv_path}[/]\n")


if __name__ == "__main__":
    main()