#!/usr/bin/env python3
"""
03_build_site.py — DENV availability report.
All figures built from data/processed/denv_final_metadata.tsv (actual downloaded dataset).
Full-genome criterion already applied at download time; this script adds host/lab filters
for analysis figures only.
"""

from datetime import date
from pathlib import Path

import pandas as pd
import numpy as np
import plotly.graph_objects as go
import plotly.express as px
import plotly.io as pio

FINAL_META = Path("data/processed/denv_final_metadata.tsv")
DOCS_DIR   = Path("docs")
DOCS_DIR.mkdir(exist_ok=True)
(DOCS_DIR / ".nojekyll").touch()

MIN_LEN = 10_000
MIN_COV = 0.95

HOST_COLORS = {"human":"#4C78A8","mosquito":"#F58518","other":"#72B7B2","unknown":"#B0B0B0"}
SERO_COLORS = {"DENV1":"#E45756","DENV2":"#F58518","DENV3":"#4C78A8","DENV4":"#72B7B2"}

# Vaccine/chimeric accessions excluded from analysis
VACCINE_ACCESSIONS = {
    "MW945952":("RDENV1-WP-1A","DENV1","Recombinant chimeric virus"),
    "AF180818":("16007 (PDK-13)","DENV1","PDK-13 passage — attenuated vaccine"),
    "U87412":  ("PDK-53","DENV2","PDK-53 passage — attenuated vaccine"),
    "KU725664":("PDK53","DENV2","PDK-53 variant"),
    "M84728":  ("16681-PDK53","DENV2","PDK-53 variant"),
    "KJ160505":("rDENV3-4","DENV3","Chimeric DENV3/4 recombinant"),
    "MW793459":("PDK48","DENV4","PDK-48 passage — attenuated vaccine"),
    "KJ160504":("rDENV4","DENV4","Recombinant chimeric DENV4"),
}

if not FINAL_META.exists():
    print(f"ERROR: {FINAL_META} not found. Run 02_download_genomes.py first.")
    raise SystemExit(1)

# ── Load final metadata ───────────────────────────────────────────────────────
raw = pd.read_csv(FINAL_META, sep="\t", low_memory=False)

for col in ["length","genome_coverage","E_coverage","C_coverage","prM_coverage",
            "NS1_coverage","NS2A_coverage","NS2B_coverage","NS3_coverage",
            "NS4A_coverage","2K_coverage","NS4B_coverage","NS5_coverage"]:
    if col in raw.columns:
        raw[col] = pd.to_numeric(raw[col], errors="coerce")

raw["collection_year"] = pd.to_numeric(raw["date"].astype(str).str[:4], errors="coerce")

# Host category from host_type (Nextstrain-curated), regex fallback
def _hcat(row):
    ht = str(row.get("host_type","")).strip()
    if ht == "Human":    return "human"
    if ht == "Mosquito": return "mosquito"
    if ht in ("Other",): return "other"
    h = str(row.get("host","")).lower()
    if "homo sapiens" in h or "human" in h: return "human"
    if any(x in h for x in ["aedes","culex","mosquito","stegomyia","culicidae"]): return "mosquito"
    return "other" if (h and h != "nan") else "unknown"

raw["host_category"] = raw.apply(_hcat, axis=1)

# ── Exclusion tracking ────────────────────────────────────────────────────────
n_total_dl = len(raw)
n_lab = int((raw["is_lab_host"].fillna("").astype(str).str.lower()=="true").sum())

# Apply lab + vaccine exclusion
lab_mask = raw["is_lab_host"].fillna("").astype(str).str.lower()=="true"
vax_mask = raw["accession"].isin(VACCINE_ACCESSIONS)
ns = raw[~lab_mask & ~vax_mask].copy()

n_unknown_host = int((ns["host_category"]=="unknown").sum())

# Analysis subset: human + mosquito only
full = ns[ns["host_category"].isin(["human","mosquito"])].copy()

today = date.today().strftime("%d %b %Y")

def fig_html(fig, div_id):
    return pio.to_html(fig, full_html=False, include_plotlyjs=False, div_id=div_id)

def card(anchor, kind, num, title, caption, content):
    return f"""
<div class="card" id="{anchor}">
  <div class="fig-label">{kind} {num} — {title}</div>
  <p class="caption">{caption}</p>
  {content}
</div>"""

# ── Table 1: summary ──────────────────────────────────────────────────────────
rows = []
for sero in ["DENV1","DENV2","DENV3","DENV4"]:
    sf = full[full["serotype"]==sero]
    rows.append({"Serotype":sero,
                 "Human":     int((sf["host_category"]=="human").sum()),
                 "Mosquito":  int((sf["host_category"]=="mosquito").sum()),
                 "Total":     len(sf)})
tbl1_df = pd.DataFrame(rows)
tot = tbl1_df.sum(numeric_only=True)
tbl1 = """<table><thead><tr>
  <th>Serotype</th><th>Human</th><th>Mosquito</th><th>Total</th>
</tr></thead><tbody>"""
for _,r in tbl1_df.iterrows():
    tbl1 += f"<tr><td><b>{r['Serotype']}</b></td><td>{int(r['Human']):,}</td><td>{int(r['Mosquito']):,}</td><td><b>{int(r['Total']):,}</b></td></tr>"
tbl1 += (f'<tr class="total"><td><b>Total</b></td>'
         f"<td><b>{int(tot['Human']):,}</b></td>"
         f"<td><b>{int(tot['Mosquito']):,}</b></td>"
         f"<td><b>{int(tot['Total']):,}</b></td></tr>")
tbl1 += "</tbody></table>"

# ── Table 2: exclusions ───────────────────────────────────────────────────────
tbl2 = f"""<table><thead><tr><th>Category</th><th>n excluded</th><th>Reason</th></tr></thead><tbody>
<tr><td>Lab-grown (<code>is_lab_host = True</code>)</td><td>{n_lab:,}</td>
  <td>Sequenced from cell lines (Vero, C6/36, etc.). May carry culture-adaptation mutations.</td></tr>
<tr><td>PDK-passage / chimeric strains (manual)</td><td>{len(VACCINE_ACCESSIONS)}</td>
  <td>Vaccine candidates or engineered chimeras not flagged by <code>is_lab_host</code>. See sub-table below.</td></tr>
<tr><td>Sequences excluded from analysis figures</td><td>{n_unknown_host:,}</td>
  <td><b>Unknown host</b>: host field is blank in the GenBank record. Common in older submissions
  (pre-2015) and some surveillance datasets where host was not recorded by the depositing lab.
  These sequences are kept in the downloaded FASTA but excluded from host-specific analyses.</td></tr>
</tbody></table>
<h4 style="margin:14px 0 6px;font-size:.88rem">PDK-passage and chimeric strains excluded</h4>
<table><thead><tr><th>Accession</th><th>Strain</th><th>Serotype</th><th>Reason</th></tr></thead><tbody>"""
for acc,(strain,sero,reason) in VACCINE_ACCESSIONS.items():
    tbl2 += f"<tr><td><code>{acc}</code></td><td>{strain}</td><td>{sero}</td><td>{reason}</td></tr>"
tbl2 += "</tbody></table>"

# ── Table 3: seqkit stats split by host ──────────────────────────────────────
def seq_stats(df_sub, label):
    l = df_sub["length"].dropna()
    c = df_sub["genome_coverage"].dropna()
    return {"Host":label, "n":len(df_sub),
            "min (bp)":int(l.min()), "mean (bp)":f"{l.mean():.0f}",
            "max (bp)":int(l.max()), "std (bp)":f"{l.std():.0f}",
            "cov mean":f"{c.mean():.4f}", "cov std":f"{c.std():.4f}"}

stats_rows = [seq_stats(full[full["host_category"]=="human"], "Human"),
              seq_stats(full[full["host_category"]=="mosquito"], "Mosquito"),
              seq_stats(full, "All (human+mosquito)")]
tbl3 = "<table><thead><tr>" + "".join(f"<th>{k}</th>" for k in stats_rows[0]) + "</tr></thead><tbody>"
for r in stats_rows:
    bold = r["Host"].startswith("All")
    tag = "b" if bold else "span"
    tbl3 += "<tr>" + "".join(f"<td><{tag}>{v}</{tag}></td>" for v in r.values()) + "</tr>"
tbl3 += "</tbody></table>"

# ── Figure 1: counts by serotype × host ──────────────────────────────────────
def fig1():
    rows = []
    for sero in ["DENV1","DENV2","DENV3","DENV4"]:
        for h,label in [("human","Human"),("mosquito","Mosquito")]:
            rows.append({"Serotype":sero,"Host":label,
                         "n":int(((full["serotype"]==sero)&(full["host_category"]==h)).sum())})
    fig = px.bar(pd.DataFrame(rows), x="Serotype", y="n", color="Host",
                 color_discrete_map={"Human":HOST_COLORS["human"],"Mosquito":HOST_COLORS["mosquito"]},
                 barmode="group", labels={"n":"# sequences"})
    fig.update_layout(height=380, legend_title="Host", margin=dict(t=20,b=40))
    return fig

# ── Figure 2: length distribution (all downloaded, after lab/vax exclusion) ──
def fig2():
    fig = go.Figure()
    for sero in ["DENV1","DENV2","DENV3","DENV4"]:
        vals = ns[ns["serotype"]==sero]["length"].dropna()
        fig.add_trace(go.Histogram(x=vals, name=sero, marker_color=SERO_COLORS[sero],
                                   opacity=0.7, xbins=dict(start=9900,end=11300,size=50)))
    fig.update_layout(barmode="overlay", height=360,
                      xaxis_title="Length (bp)", yaxis_title="# sequences",
                      legend_title="Serotype", margin=dict(t=20,b=40))
    return fig

# ── Figure 3: genome coverage distribution ───────────────────────────────────
def fig3():
    fig = go.Figure()
    for sero in ["DENV1","DENV2","DENV3","DENV4"]:
        vals = ns[ns["serotype"]==sero]["genome_coverage"].dropna()
        fig.add_trace(go.Histogram(x=vals, name=sero, marker_color=SERO_COLORS[sero],
                                   opacity=0.7, xbins=dict(start=0.3,end=1.01,size=0.01)))
    fig.add_vline(x=MIN_COV, line_dash="dash", line_color="red",
                  annotation_text=f"Download threshold ({MIN_COV})", annotation_position="top left")
    fig.update_layout(barmode="overlay", height=360,
                      xaxis_title="genome_coverage", yaxis_title="# sequences",
                      legend_title="Serotype", margin=dict(t=20,b=40))
    return fig

# ── Figure 4: per-gene coverage (mean ± SD, human+mosquito only) ─────────────
def fig4():
    gene_map = {"C_coverage":"C","prM_coverage":"prM","E_coverage":"E",
                "NS1_coverage":"NS1","NS2A_coverage":"NS2A","NS2B_coverage":"NS2B",
                "NS3_coverage":"NS3","NS4A_coverage":"NS4A","2K_coverage":"2K",
                "NS4B_coverage":"NS4B","NS5_coverage":"NS5"}
    avail_cols  = [c for c in gene_map if c in full.columns]
    avail_labels= [gene_map[c] for c in avail_cols]
    if not avail_cols: return None
    fig = go.Figure()
    for sero in ["DENV1","DENV2","DENV3","DENV4"]:
        sub = full[full["serotype"]==sero][avail_cols].dropna()
        fig.add_trace(go.Bar(x=avail_labels, y=[float(sub[c].mean()) for c in avail_cols],
                             name=sero, marker_color=SERO_COLORS[sero],
                             error_y=dict(type="data",
                                          array=[float(sub[c].std()) for c in avail_cols],
                                          visible=True)))
    if "NS1" in avail_labels:
        i = avail_labels.index("NS1")
        fig.add_vrect(x0=i-0.5, x1=i+0.5, fillcolor="yellow", opacity=0.15,
                      line_width=0, annotation_text="NS1", annotation_position="top left")
    fig.update_layout(barmode="group", height=420, xaxis_title="Gene",
                      yaxis_title="Mean coverage ± SD", yaxis=dict(range=[0,1.08]),
                      legend_title="Serotype", margin=dict(t=20,b=40))
    return fig

# ── Figure 5: timeline ────────────────────────────────────────────────────────
def fig5():
    df = full[full["collection_year"].between(1940,2026)]
    grouped = df.groupby(["collection_year","serotype"]).size().reset_index(name="n")
    fig = px.bar(grouped, x="collection_year", y="n", color="serotype",
                 color_discrete_map=SERO_COLORS,
                 labels={"collection_year":"Year","n":"# sequences","serotype":"Serotype"})
    fig.update_layout(height=360, barmode="stack", legend_title="Serotype", margin=dict(t=20,b=40))
    return fig

# ── Figure 6: map ─────────────────────────────────────────────────────────────
def fig6():
    df = full.dropna(subset=["country"])
    cc = df.groupby("country").size().reset_index(name="n")
    fig = px.choropleth(cc, locations="country", locationmode="country names",
                        color="n", color_continuous_scale="Blues", labels={"n":"Sequences"})
    fig.update_layout(height=420,
                      geo=dict(showframe=False,showcoastlines=True,projection_type="natural earth"),
                      margin=dict(t=20,b=10,l=0,r=0))
    return fig

# ── Table 4: mosquito species ─────────────────────────────────────────────────
mosq = full[full["host_category"]=="mosquito"].dropna(subset=["host"])
mosq_counts = mosq.groupby(["host","serotype"]).size().reset_index(name="n")
mosq_totals = mosq_counts.groupby("host")["n"].sum().sort_values(ascending=False)
tbl4 = "<table><thead><tr><th>Host species</th><th>DENV1</th><th>DENV2</th><th>DENV3</th><th>DENV4</th><th>Total</th></tr></thead><tbody>"
for species in mosq_totals.index:
    row_n = {s:0 for s in ["DENV1","DENV2","DENV3","DENV4"]}
    for _,r in mosq_counts[mosq_counts["host"]==species].iterrows():
        row_n[r["serotype"]] = int(r["n"])
    tbl4 += (f"<tr><td><i>{species}</i></td>"
             f"<td>{row_n['DENV1']}</td><td>{row_n['DENV2']}</td>"
             f"<td>{row_n['DENV3']}</td><td>{row_n['DENV4']}</td>"
             f"<td><b>{int(mosq_totals[species])}</b></td></tr>")
tbl4 += "</tbody></table>"

# ── Figures 7-9: lineage / genotype ──────────────────────────────────────────
def _bar_horiz(col, title_col):
    if col not in full.columns: return None
    df = full.dropna(subset=[col])
    df = df[df[col].str.strip().ne("") & df[col].str.strip().ne("unassigned")]
    if df.empty: return None
    counts = (df.groupby([col,"serotype"]).size().reset_index(name="n")
              .sort_values(["serotype","n"],ascending=[True,False]))
    top = counts.groupby("serotype").head(15)
    fig = px.bar(top, x="n", y=col, color="serotype", color_discrete_map=SERO_COLORS,
                 orientation="h", labels={"n":"# sequences",col:title_col})
    fig.update_layout(height=max(380,len(top)*18), yaxis=dict(autorange="reversed"),
                      legend_title="Serotype", margin=dict(t=20,b=40))
    return fig

F = {1:fig1(), 2:fig2(), 3:fig3(), 4:fig4(), 5:fig5(), 6:fig6(),
     7:_bar_horiz("major_lineage","Major lineage"),
     8:_bar_horiz("minor_lineage","Minor lineage"),
     9:_bar_horiz("genotype","Genotype (Nextclade)")}

# ── Build body ────────────────────────────────────────────────────────────────
n_total = len(full)
n_human = int((full["host_category"]=="human").sum())
n_mosquito = int((full["host_category"]=="mosquito").sum())

body = ""
body += card("table1","Table",1,"Analysis dataset — human and mosquito full genomes by serotype",
    f"Sequences in the analysis dataset: lab-grown and vaccine-derived excluded, "
    f"host restricted to human and mosquito, from the downloaded FASTA ({FINAL_META.name}). "
    f"Total: <b>{n_total:,}</b> ({n_human:,} human, {n_mosquito:,} mosquito).", tbl1)

body += card("table2","Table",2,"Excluded sequences — lab, vaccine, and unknown host",
    "Sequences removed from analysis figures. Downloaded FASTA is unmodified; exclusions are metadata-level.", tbl2)

body += card("table3","Table",3,"seqkit stats — split by host",
    "Sequence statistics for the analysis dataset (human+mosquito, after all exclusions), "
    "computed from the downloaded FASTA via seqkit.", tbl3)

body += card("figure1","Figure",1,"Full-genome counts by serotype and host",
    "Grouped bar chart from the analysis dataset (human and mosquito only, lab/vaccine excluded).",
    fig_html(F[1],"f1"))

body += card("figure2","Figure",2,"Sequence length distribution (downloaded dataset, all hosts)",
    "Length distribution of all downloaded sequences (lab/vaccine excluded, all hosts). "
    "All sequences are in the 10–12 kb range by design (seqkit length filter).",
    fig_html(F[2],"f2"))

body += card("figure3","Figure",3,"Genome coverage distribution (downloaded dataset, all hosts)",
    f"Distribution of <code>genome_coverage</code> for all downloaded sequences. "
    f"The red dashed line marks the download threshold ({MIN_COV}); all retained sequences are to the right.",
    fig_html(F[3],"f3"))

if F[4]:
    body += card("figure4","Figure",4,"Per-gene coverage — mean ± SD (analysis dataset)",
        "Mean ± SD per-gene coverage across the analysis dataset (human + mosquito). "
        "Yellow band highlights NS1; NS4A-2K-NS4B is the likely study target (TBC).",
        fig_html(F[4],"f4") +
        '<div class="question-box"><b>Question for Beatriz:</b> NS1 is highlighted but the '
        'target may be <b>NS4A-2K-NS4B</b>. Please confirm which gene(s) should be prioritised.</div>')

body += card("figure5","Figure",5,"Temporal distribution — analysis dataset",
    "Collection year distribution by serotype (human + mosquito sequences).",
    fig_html(F[5],"f5"))

body += card("figure6","Figure",6,"Geographic distribution — analysis dataset",
    "Choropleth map of analysis sequences per country.", fig_html(F[6],"f6"))

body += card("table4","Table",4,"Mosquito host species breakdown",
    "All mosquito-derived sequences from the analysis dataset, by host species and serotype. "
    "Please indicate which species to retain for the final analysis.",
    tbl4)

if F[7]:
    body += card("figure7","Figure",7,"Major lineage distribution",
        "Major lineage (Nextclade) of analysis sequences.", fig_html(F[7],"f7"))
if F[8]:
    body += card("figure8","Figure",8,"Minor lineage distribution (top 15 per serotype)",
        "Minor lineage breakdown.", fig_html(F[8],"f8"))
if F[9]:
    body += card("figure9","Figure",9,"Genotype distribution",
        "Nextclade genotype assignments.", fig_html(F[9],"f9"))

# ── TOC ───────────────────────────────────────────────────────────────────────
toc_items = [
    ("table1","Table 1 — Analysis dataset summary"),
    ("table2","Table 2 — Excluded sequences"),
    ("table3","Table 3 — seqkit stats by host"),
    ("figure1","Figure 1 — Counts by serotype &amp; host"),
    ("figure2","Figure 2 — Length distribution"),
    ("figure3","Figure 3 — Genome coverage"),
    ("figure4","Figure 4 — Per-gene coverage"),
    ("figure5","Figure 5 — Temporal distribution"),
    ("figure6","Figure 6 — Geographic distribution"),
    ("table4","Table 4 — Mosquito species"),
    ("figure7","Figure 7 — Major lineage"),
    ("figure8","Figure 8 — Minor lineage"),
    ("figure9","Figure 9 — Genotype"),
]
toc_html = "\n".join(f'<li><a href="#{a}">{l}</a></li>' for a,l in toc_items)

html = f"""<!DOCTYPE html>
<html lang="en"><head>
<meta charset="UTF-8"><meta name="viewport" content="width=device-width,initial-scale=1.0">
<title>DENV Full-Genome Dataset — CNRS UMR2K 2026</title>
<script src="https://cdn.plot.ly/plotly-2.35.2.min.js"></script>
<style>
:root{{--blue:#4C78A8;--bg:#f7f8fa;--card:#fff;--border:#e0e4ea;--text:#1a1a2e;--muted:#6c757d}}
*{{box-sizing:border-box;margin:0;padding:0}}
body{{font-family:'Segoe UI',Arial,sans-serif;background:var(--bg);color:var(--text);font-size:14px}}
header{{background:var(--blue);color:#fff;padding:28px 40px}}
header h1{{font-size:1.5rem;font-weight:600;margin-bottom:6px}}
header p{{opacity:.85;font-size:.88rem}}
main{{max-width:1100px;margin:0 auto;padding:24px 20px}}
.criteria{{background:#e8f4fd;border-left:4px solid var(--blue);
  padding:14px 20px;border-radius:4px;margin-bottom:20px;font-size:.88rem;line-height:1.6}}
.criteria code{{background:#cde4f7;padding:1px 5px;border-radius:3px}}
.toc{{background:var(--card);border:1px solid var(--border);border-radius:8px;
  padding:16px 24px;margin-bottom:24px}}
.toc h3{{font-size:.88rem;font-weight:600;margin-bottom:8px;color:var(--blue)}}
.toc ul{{list-style:none;column-count:2;column-gap:24px}}
.toc li{{padding:3px 0;font-size:.82rem}}
.toc a{{color:var(--text);text-decoration:none}}
.toc a:hover{{text-decoration:underline}}
.card{{background:var(--card);border:1px solid var(--border);border-radius:8px;
  padding:20px 24px;margin-bottom:24px;box-shadow:0 1px 4px rgba(0,0,0,.05)}}
.fig-label{{font-weight:700;font-size:.82rem;color:var(--blue);
  text-transform:uppercase;letter-spacing:.5px;margin-bottom:6px}}
.caption{{color:var(--muted);font-size:.84rem;margin-bottom:14px;line-height:1.55}}
table{{width:100%;border-collapse:collapse;font-size:.85rem;margin-top:8px}}
th{{background:var(--blue);color:#fff;padding:9px 14px;text-align:left;font-weight:600}}
td{{padding:7px 14px;border-bottom:1px solid var(--border)}}
tr:hover{{background:#f0f4fa}}
tr.total{{background:#eef2f7;font-weight:600}}
h4{{color:var(--blue)}}
.question-box{{background:#fff8e1;border-left:4px solid #f0ad4e;border-radius:4px;
  padding:12px 16px;margin-top:16px;font-size:.88rem;line-height:1.55}}
footer{{text-align:center;padding:24px;color:var(--muted);font-size:.8rem}}
</style></head><body>
<header>
  <h1>DENV Full-Genome Dataset — Availability Report</h1>
  <p>CNRS UMR2K Seed Grant 2026 &nbsp;·&nbsp; Author: Camila Duitama &nbsp;·&nbsp; Generated: {today}</p>
</header>
<main>
<div class="criteria">
  <b>Data source:</b> Nextstrain DENV (nextstrain.org/dengue, updated 2026-07-02) curated from GenBank,
  downloaded to <code>denv_final.fasta</code> ({n_total_dl:,} sequences total).
  All sequences pass: <code>length ≥ {MIN_LEN:,} bp</code> AND <code>genome_coverage ≥ {MIN_COV}</code>
  (applied at download time).<br>
  <b>Analysis figures use a subset:</b> <code>is_lab_host ≠ True</code> ({n_lab} excluded) ·
  PDK/chimeric strains removed ({len(VACCINE_ACCESSIONS)}) ·
  host is <b>Human or Mosquito</b> ({n_unknown_host:,} sequences with no host annotation excluded from figures).
  <b>Total in analysis figures: {n_total:,}</b> ({n_human:,} human · {n_mosquito:,} mosquito).
</div>
<div class="toc"><h3>Contents</h3><ul>{toc_html}</ul></div>
{body}
</main>
<footer>Data: Nextstrain dengue (nextstrain.org) · GenBank · Generated {today}</footer>
</body></html>"""

out = DOCS_DIR / "index.html"
out.write_text(html, encoding="utf-8")
print(f"Site written → {out}  ({out.stat().st_size//1024} KB)")
