#!/usr/bin/env python3
"""
03_build_site.py — DENV availability report (Nextstrain-only, Nextstrain-first).
Full-genome criterion: length >= 10,000 bp AND genome_coverage >= 0.95 AND is_lab_host != True
"""

from datetime import date
from pathlib import Path

import pandas as pd
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
import plotly.io as pio

RAW_DIR  = Path("data/raw")
DOCS_DIR = Path("docs")
DOCS_DIR.mkdir(exist_ok=True)
(DOCS_DIR / ".nojekyll").touch()

MIN_LEN = 10_000
MIN_COV = 0.95

HOST_COLORS  = {"human":"#4C78A8","mosquito":"#F58518","other":"#72B7B2","unknown":"#B0B0B0"}
SERO_COLORS  = {"DENV1":"#E45756","DENV2":"#F58518","DENV3":"#4C78A8","DENV4":"#72B7B2"}

# ── Load all Nextstrain metadata ──────────────────────────────────────────────
dfs = []
for sero in ["denv1","denv2","denv3","denv4"]:
    p = RAW_DIR / f"nextstrain_metadata_{sero}.tsv"
    if not p.exists():
        print(f"WARNING: {p} missing — run 01_query_availability.py first")
        continue
    df = pd.read_csv(p, sep="\t", low_memory=False)
    df["serotype"] = sero.upper()
    dfs.append(df)

ns = pd.concat(dfs, ignore_index=True)

for col in ["length","genome_coverage","NS1_coverage","E_coverage","prM_coverage",
            "NS2A_coverage","NS2B_coverage","NS3_coverage","NS4A_coverage",
            "NS4B_coverage","NS5_coverage","2K_coverage"]:
    if col in ns.columns:
        ns[col] = pd.to_numeric(ns[col], errors="coerce")

ns["collection_year"] = pd.to_numeric(ns["date"].astype(str).str[:4], errors="coerce")

# Host category — use host_type (Nextstrain-curated) when available, regex fallback
def _host_cat(row):
    ht = str(row.get("host_type","")).strip()
    if ht == "Human":   return "human"
    if ht == "Mosquito": return "mosquito"
    if ht in ("Other","Mus"): return "other"
    # regex fallback
    h = str(row.get("host","")).lower()
    if "homo sapiens" in h or "human" in h: return "human"
    if any(x in h for x in ["aedes","culex","mosquito","stegomyia","culicidae"]): return "mosquito"
    if h and h != "nan": return "other"
    return "unknown"

ns["host_category"] = ns.apply(_host_cat, axis=1)

# Exclude lab sequences
lab = ns["is_lab_host"].fillna("").astype(str).str.lower() == "true"
n_lab_excluded = int(lab.sum())
ns_nolab = ns[~lab].copy()

# PDK/chimeric strains excluded by name (not caught by is_lab_host)
VACCINE_ACCESSIONS = {
    "MW945952": ("RDENV1-WP-1A",    "DENV1", "Recombinant chimeric virus"),
    "AF180818": ("16007 (PDK-13)",  "DENV1", "PDK-13 passage — attenuated vaccine candidate"),
    "U87412":   ("PDK-53",          "DENV2", "PDK-53 passage — attenuated vaccine candidate"),
    "KU725664": ("PDK53",           "DENV2", "PDK-53 variant"),
    "M84728":   ("16681-PDK53",     "DENV2", "PDK-53 variant"),
    "KJ160505": ("rDENV3-4",        "DENV3", "Chimeric DENV3/4 recombinant"),
    "MW793459": ("PDK48",           "DENV4", "PDK-48 passage — attenuated vaccine candidate"),
    "KJ160504": ("rDENV4",          "DENV4", "Recombinant chimeric DENV4"),
}
vax_mask = ns_nolab["accession"].isin(VACCINE_ACCESSIONS)
ns_nolab = ns_nolab[~vax_mask].copy()

# Full-genome filter — human and mosquito only
full = ns_nolab[
    (ns_nolab["length"] >= MIN_LEN) &
    (ns_nolab["genome_coverage"] >= MIN_COV) &
    (ns_nolab["host_category"].isin(["human", "mosquito"]))
].copy()

def fig_html(fig, div_id):
    return pio.to_html(fig, full_html=False, include_plotlyjs=False, div_id=div_id)

# ── Table 1: Summary counts ───────────────────────────────────────────────────
rows = []
for sero in ["DENV1","DENV2","DENV3","DENV4"]:
    sf = full[full["serotype"]==sero]
    rows.append({
        "Serotype":  sero,
        "Human":     int((sf["host_category"]=="human").sum()),
        "Mosquito":  int((sf["host_category"]=="mosquito").sum()),
        "Total":     len(sf),
    })
tbl1_df = pd.DataFrame(rows)

tbl1 = """<table>
<thead><tr>
  <th>Serotype</th>
  <th>Human</th><th>Mosquito</th><th>Total</th>
</tr></thead><tbody>
"""
for _, r in tbl1_df.iterrows():
    tbl1 += (
        f"<tr><td><b>{r['Serotype']}</b></td>"
        f"<td>{int(r['Human']):,}</td>"
        f"<td>{int(r['Mosquito']):,}</td>"
        f"<td><b>{int(r['Total']):,}</b></td></tr>\n"
    )
tot = tbl1_df.sum(numeric_only=True)
tbl1 += (
    f'<tr class="total"><td><b>Total</b></td>'
    f"<td><b>{int(tot['Human']):,}</b></td>"
    f"<td><b>{int(tot['Mosquito']):,}</b></td>"
    f"<td><b>{int(tot['Total']):,}</b></td></tr>"
)
tbl1 += "</tbody></table>"

# ── Table 2: Mosquito species breakdown ──────────────────────────────────────
mosq_full = full[full["host_category"]=="mosquito"].dropna(subset=["host"])
mosq_counts = mosq_full.groupby(["host","serotype"]).size().reset_index(name="n")
mosq_totals = mosq_counts.groupby("host")["n"].sum().sort_values(ascending=False)

tbl2 = """<table>
<thead><tr>
  <th>Host species</th>
  <th>DENV1</th><th>DENV2</th><th>DENV3</th><th>DENV4</th><th>Total</th>
</tr></thead><tbody>
"""
for species in mosq_totals.index:
    row_n = {s:0 for s in ["DENV1","DENV2","DENV3","DENV4"]}
    sub = mosq_counts[mosq_counts["host"]==species]
    for _, r in sub.iterrows():
        row_n[r["serotype"]] = int(r["n"])
    tbl2 += (
        f"<tr><td><i>{species}</i></td>"
        f"<td>{row_n['DENV1']}</td><td>{row_n['DENV2']}</td>"
        f"<td>{row_n['DENV3']}</td><td>{row_n['DENV4']}</td>"
        f"<td><b>{int(mosq_totals[species])}</b></td></tr>\n"
    )
tbl2 += "</tbody></table>"

# ── Figure 1: Full-genome counts by serotype × host ──────────────────────────
def fig1():
    rows = []
    for sero in ["DENV1","DENV2","DENV3","DENV4"]:
        for hcat, label in [("human","Human"),("mosquito","Mosquito")]:
            n = int(((full["serotype"]==sero)&(full["host_category"]==hcat)).sum())
            rows.append({"Serotype":sero,"Host":label,"n":n})
    df = pd.DataFrame(rows)
    fig = px.bar(df, x="Serotype", y="n", color="Host",
                 color_discrete_map={"Human":HOST_COLORS["human"],"Mosquito":HOST_COLORS["mosquito"]},
                 barmode="group",
                 labels={"n":"# full genomes","Host":"Host"})
    fig.update_layout(height=400, legend_title="Host", margin=dict(t=20,b=40))
    return fig

# ── Figure 2: Length distribution ────────────────────────────────────────────
def fig2():
    fig = go.Figure()
    for sero in ["DENV1","DENV2","DENV3","DENV4"]:
        sub = ns_nolab[ns_nolab["serotype"]==sero]["length"].dropna()
        fig.add_trace(go.Histogram(x=sub, name=sero, marker_color=SERO_COLORS[sero],
                                   opacity=0.7, xbins=dict(start=0,end=12000,size=200)))
    fig.add_vline(x=MIN_LEN, line_dash="dash", line_color="red",
                  annotation_text="10 kb length threshold", annotation_position="top right")
    fig.update_layout(barmode="overlay", height=380,
                      xaxis_title="Length (bp)", yaxis_title="# sequences",
                      legend_title="Serotype", margin=dict(t=20,b=40))
    return fig

# ── Figure 3: Genome coverage distribution (≥10kb sequences) ─────────────────
def fig3():
    sub = ns_nolab[ns_nolab["length"]>=MIN_LEN].dropna(subset=["genome_coverage"])
    fig = go.Figure()
    for sero in ["DENV1","DENV2","DENV3","DENV4"]:
        vals = sub[sub["serotype"]==sero]["genome_coverage"]
        fig.add_trace(go.Histogram(x=vals, name=sero, marker_color=SERO_COLORS[sero],
                                   opacity=0.7, xbins=dict(start=0,end=1.01,size=0.01)))
    fig.add_vline(x=MIN_COV, line_dash="dash", line_color="red",
                  annotation_text=f"Coverage threshold ({MIN_COV})",
                  annotation_position="top left")
    fig.update_layout(barmode="overlay", height=380,
                      xaxis_title="genome_coverage (fraction of reference covered)",
                      yaxis_title="# sequences (length ≥10 kb)",
                      legend_title="Serotype", margin=dict(t=20,b=40))
    return fig

# ── Figure 4: Per-gene coverage — violin/box plot ───────────────────────────
def fig4():
    gene_map = {"C_coverage":"C","prM_coverage":"prM","E_coverage":"E",
                "NS1_coverage":"NS1","NS2A_coverage":"NS2A","NS2B_coverage":"NS2B",
                "NS3_coverage":"NS3","NS4A_coverage":"NS4A","2K_coverage":"2K",
                "NS4B_coverage":"NS4B","NS5_coverage":"NS5"}
    avail_cols  = [c for c in gene_map if c in full.columns]
    avail_labels= [gene_map[c] for c in avail_cols]
    if not avail_cols: return None

    fig = go.Figure()
    import numpy as np
    for sero in ["DENV1","DENV2","DENV3","DENV4"]:
        sub = full[full["serotype"]==sero][avail_cols].dropna()
        fig.add_trace(go.Bar(
            x=avail_labels,
            y=[float(sub[c].mean()) for c in avail_cols],
            name=sero,
            marker_color=SERO_COLORS[sero],
            error_y=dict(
                type="data",
                array=[float(sub[c].std()) for c in avail_cols],
                visible=True,
            ),
        ))
    # Highlight NS1
    if "NS1" in avail_labels:
        i = avail_labels.index("NS1")
        fig.add_vrect(x0=i-0.5, x1=i+0.5, fillcolor="yellow", opacity=0.15,
                      line_width=0, annotation_text="NS1 (study target)",
                      annotation_position="top left")
    fig.update_layout(barmode="group", height=420,
                      xaxis_title="Gene", yaxis_title="Mean coverage ± SD",
                      yaxis=dict(range=[0,1.08]), legend_title="Serotype",
                      margin=dict(t=20,b=40))
    return fig

# ── Figure 5: Timeline ────────────────────────────────────────────────────────
def fig5():
    df = full[full["collection_year"].between(1940,2026)]
    grouped = df.groupby(["collection_year","serotype"]).size().reset_index(name="n")
    fig = px.bar(grouped, x="collection_year", y="n", color="serotype",
                 color_discrete_map=SERO_COLORS,
                 labels={"collection_year":"Year","n":"# sequences","serotype":"Serotype"})
    fig.update_layout(height=380, barmode="stack", legend_title="Serotype",
                      margin=dict(t=20,b=40))
    return fig

# ── Figure 6: Geographic map ──────────────────────────────────────────────────
def fig6():
    df = full.dropna(subset=["country"])
    cc = df.groupby("country").size().reset_index(name="n")
    fig = px.choropleth(cc, locations="country", locationmode="country names",
                        color="n", color_continuous_scale="Blues",
                        labels={"n":"Full genomes"})
    fig.update_layout(height=440,
                      geo=dict(showframe=False,showcoastlines=True,
                               projection_type="natural earth"),
                      margin=dict(t=20,b=10,l=0,r=0),
                      coloraxis_colorbar_title="# sequences")
    return fig

# ── Figure 7: Major lineage ───────────────────────────────────────────────────
def fig7():
    col = "major_lineage"
    if col not in full.columns: return None
    df = full.dropna(subset=[col])
    df = df[df[col].str.strip().ne("") & df[col].str.strip().ne("unassigned")]
    counts = (df.groupby([col,"serotype"]).size().reset_index(name="n")
              .sort_values(["serotype","n"],ascending=[True,False]))
    if counts.empty: return None
    fig = px.bar(counts, x="n", y=col, color="serotype",
                 color_discrete_map=SERO_COLORS, orientation="h",
                 labels={"n":"# full genomes",col:"Major lineage"})
    fig.update_layout(height=max(400,len(counts)*18),
                      yaxis=dict(autorange="reversed"),
                      legend_title="Serotype", margin=dict(t=20,b=40))
    return fig

# ── Figure 8: Minor lineage ───────────────────────────────────────────────────
def fig8():
    col = "minor_lineage"
    if col not in full.columns: return None
    df = full.dropna(subset=[col])
    df = df[df[col].str.strip().ne("") & df[col].str.strip().ne("unassigned")]
    if df.empty: return None
    counts = (df.groupby([col,"serotype"]).size().reset_index(name="n")
              .sort_values(["serotype","n"],ascending=[True,False]))
    top = counts.groupby("serotype").head(15)  # top 15 per serotype
    fig = px.bar(top, x="n", y=col, color="serotype",
                 color_discrete_map=SERO_COLORS, orientation="h",
                 labels={"n":"# full genomes",col:"Minor lineage"})
    fig.update_layout(height=max(400,len(top)*18),
                      yaxis=dict(autorange="reversed"),
                      legend_title="Serotype", margin=dict(t=20,b=40))
    return fig

# ── Figure 9: Genotype distribution ──────────────────────────────────────────
def fig9():
    col = "genotype"
    if col not in full.columns: return None
    df = full.dropna(subset=[col])
    df = df[df[col].str.strip().ne("")]
    counts = (df.groupby([col,"serotype"]).size().reset_index(name="n")
              .sort_values(["serotype","n"],ascending=[True,False]))
    if counts.empty: return None
    fig = px.bar(counts, x="n", y=col, color="serotype",
                 color_discrete_map=SERO_COLORS, orientation="h",
                 labels={"n":"# full genomes",col:"Genotype (Nextclade)"})
    fig.update_layout(height=max(400,len(counts)*22),
                      yaxis=dict(autorange="reversed"),
                      legend_title="Serotype", margin=dict(t=20,b=40))
    return fig

# ── Render all ────────────────────────────────────────────────────────────────
F = {1:fig1(), 2:fig2(), 3:fig3(), 4:fig4(), 5:fig5(),
     6:fig6(), 7:fig7(), 8:fig8(), 9:fig9()}

def card(anchor, kind, num, title, caption, content):
    return f"""
<div class="card" id="{anchor}">
  <div class="fig-label">{kind} {num} — {title}</div>
  <p class="caption">{caption}</p>
  {content}
</div>"""

today = date.today().strftime("%d %b %Y")
n_total   = len(full)
n_human   = int((full["host_category"]=="human").sum())
n_mosquito= int((full["host_category"]=="mosquito").sum())

# ── Table 3: Exclusion summary ────────────────────────────────────────────────
tbl3 = f"""<table>
<thead><tr><th>Category</th><th>n excluded</th><th>Reason</th></tr></thead><tbody>
<tr>
  <td>Lab-grown sequences (<code>is_lab_host = True</code>)</td>
  <td>{n_lab_excluded:,}</td>
  <td>Sequences passaged in cell lines (Vero, C6/36, etc.); host annotation is a cell line,
  not a natural organism. May carry culture-adaptation mutations.</td>
</tr>
<tr>
  <td>PDK-passage / chimeric strains (manual)</td>
  <td>{len(VACCINE_ACCESSIONS)}</td>
  <td>Vaccine candidates or engineered chimeras not flagged by <code>is_lab_host</code>.
  Identified by scanning strain names. PDK = serially passaged in Primary Dog Kidney cells
  (attenuated). rDEN = recombinant chimeric backbone. All carry non-natural mutation loads.</td>
</tr>
<tr>
  <td>Other host (primates, rodents, unknown)</td>
  <td>—</td>
  <td>Not counted here; simply outside the human/mosquito scope of this project.</td>
</tr>
<tr class="total">
  <td><b>Total excluded (lab/vaccine)</b></td>
  <td><b>{n_lab_excluded + len(VACCINE_ACCESSIONS):,}</b></td>
  <td></td>
</tr>
</tbody></table>

<h4 style="margin:14px 0 6px;font-size:.88rem">PDK-passage and chimeric strains excluded by accession</h4>
<table>
<thead><tr><th>Accession</th><th>Strain name</th><th>Serotype</th><th>Reason</th></tr></thead>
<tbody>"""
for acc, (strain, sero, reason) in VACCINE_ACCESSIONS.items():
    tbl3 += f"<tr><td><code>{acc}</code></td><td>{strain}</td><td>{sero}</td><td>{reason}</td></tr>\n"
tbl3 += "</tbody></table>"
n_lab     = int(lab.sum())

body = ""

body += card("table1","Table",1,"Dataset summary — human and mosquito full genomes by serotype",
    f"Counts of full genomes (length ≥{MIN_LEN:,} bp and genome_coverage ≥{MIN_COV}) "
    f"from human and mosquito hosts only. Lab-adapted and vaccine-derived sequences excluded (see Table 2).",
    tbl1)

body += card("table3","Table",2,"Excluded sequences — lab-grown, vaccine-derived, and chimeric strains",
    "Sequences removed before analysis. "
    f"<code>is_lab_host = True</code> sequences are flagged by Nextstrain from their GenBank host annotation. "
    "The additional PDK-passage and chimeric strains were identified manually by scanning strain names "
    "and are not caught by the automated flag.",
    tbl3)

body += card("figure1","Figure",1,"Full-genome counts by serotype and host",
    "Grouped bar chart showing full genomes (human and mosquito only) per serotype. "
    "All other host categories are excluded from this report.",
    fig_html(F[1],"f1"))

body += card("figure2","Figure",2,"Sequence length distribution — all Nextstrain sequences",
    f"Histogram of all sequences before quality filtering. The red dashed line marks the "
    f"{MIN_LEN:,} bp length threshold. "
    "Note the bimodal distribution: short amplicons cluster below 2 kb; "
    "near-complete and complete genomes cluster above 9 kb. "
    "Sequences to the right of the cutoff are further filtered by genome coverage (Figure 3).",
    fig_html(F[2],"f2"))

body += card("figure3","Figure",3,"Genome coverage distribution — sequences ≥10 kb (all hosts)",
    f"Distribution of <code>genome_coverage</code> (fraction of the reference genome covered by the alignment) "
    f"for <b>all</b> Nextstrain sequences passing the length filter, regardless of host. "
    f"Shown here for context: the red dashed line at {MIN_COV} is the threshold applied in this report. "
    "Length alone is insufficient: a sequence can be ≥10 kb yet have low coverage due to amplicon tiling "
    "or assembly gaps. Only sequences to the <b>right</b> of the red line (and with human or mosquito host) "
    "are included in all other figures and tables.",
    fig_html(F[3],"f3"))

if F[4]:
    body += card("figure4","Figure",4,"Mean per-gene coverage across full genomes (mean ± SD)",
        "Mean coverage per gene across all human and mosquito full genomes, with error bars showing ±1 SD. "
        "Values close to 1.0 indicate that the gene is fully covered. "
        "The yellow band highlights <b>NS1</b> — the primary study target gene. "
        "Low mean or high SD for a gene may indicate incomplete assemblies in a subset of sequences.",
        fig_html(F[4],"f4") +
        '<div class="question-box">'
        '<b>Question for Beatriz:</b> NS1 is the primary gene of interest for this project. '
        'Are there <b>other genes</b> (e.g. E, prM, NS3, NS5) that should also be prioritised '
        'or verified for completeness in the final download set?'
        '</div>')

body += card("figure5","Figure",5,"Temporal distribution of full genomes",
    "Collection year distribution of full genomes by serotype. "
    "Sequences with missing or unparseable collection dates are excluded. "
    "This shows the historical sampling depth available for phylogenetic analyses.",
    fig_html(F[5],"f5"))

body += card("figure6","Figure",6,"Geographic distribution of full genomes",
    "Choropleth map coloured by number of full genomes per country. "
    "White countries have no sequences. "
    "Geographic representation directly affects the phylogeographic scope of WP1/WP2.",
    fig_html(F[6],"f6"))

body += card("table2","Table",2,"Mosquito host species — full genomes by species and serotype",
    "All mosquito-derived full genomes broken down by host species (from Nextstrain <code>host</code> field) "
    "and serotype. <i>Culicidae</i> entries are identified to family level only — species unknown. "
    "African sylvatic species (<i>A. luteocephalus</i>, <i>A. taylori</i>, <i>A. africanus</i>) "
    "are part of the ancestral sylvatic dengue cycle. "
    "Please indicate which species should be included or excluded in the final download.",
    tbl2)

if F[7]:
    body += card("figure7","Figure",7,"Major lineage distribution of full genomes",
        "Full genomes grouped by major lineage (Nextclade-assigned). "
        "Lineage diversity determines how representative the dataset is of known viral diversity. "
        "'unassigned' sequences are excluded from this figure but remain in the download set.",
        fig_html(F[7],"f7"))

if F[8]:
    body += card("figure8","Figure",8,"Minor lineage distribution of full genomes (top 15 per serotype)",
        "Top 15 minor lineages per serotype. Fine-grained lineage structure relevant to WP1/WP2 "
        "phylogenetic analyses.",
        fig_html(F[8],"f8"))

if F[9]:
    body += card("figure9","Figure",9,"Genotype distribution of full genomes",
        "Genotype assignments from Nextclade. Each DENV serotype comprises multiple genotypes with "
        "distinct geographic and evolutionary histories. This distribution informs sampling strategy "
        "for phylogenetic modelling.",
        fig_html(F[9],"f9"))

# TOC
toc_items = [
    ("table1","Table 1 — Dataset summary"),
    ("table3","Table 2 — Excluded sequences"),
    ("figure1","Figure 1 — Counts by serotype &amp; host"),
    ("figure2","Figure 2 — Length distribution"),
    ("figure3","Figure 3 — Genome coverage distribution"),
    ("figure4","Figure 4 — Per-gene coverage"),
    ("figure5","Figure 5 — Temporal distribution"),
    ("figure6","Figure 6 — Geographic distribution"),
    ("table2","Table 2 — Mosquito species &amp; relevance"),
    ("figure7","Figure 7 — Major lineage"),
    ("figure8","Figure 8 — Minor lineage"),
    ("figure9","Figure 9 — Genotype distribution"),
]
toc_html = "\n".join(f'<li><a href="#{a}">{l}</a></li>' for a,l in toc_items)

html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width,initial-scale=1.0">
<title>DENV Genome Availability — CNRS UMR2K 2026</title>
<script src="https://cdn.plot.ly/plotly-2.35.2.min.js"></script>
<style>
:root{{--blue:#4C78A8;--orange:#F58518;--teal:#72B7B2;--red:#E45756;
  --bg:#f7f8fa;--card:#ffffff;--border:#e0e4ea;--text:#1a1a2e;--muted:#6c757d}}
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
.tbl-note{{font-size:.8rem;color:var(--muted);margin-top:10px;line-height:1.5}}
.question-box{{background:#fff8e1;border-left:4px solid #f0ad4e;border-radius:4px;
  padding:12px 16px;margin-top:16px;font-size:.88rem;line-height:1.55}}
footer{{text-align:center;padding:24px;color:var(--muted);font-size:.8rem}}
</style>
</head>
<body>
<header>
  <h1>DENV Full-Genome Dataset — Availability Report</h1>
  <p>CNRS UMR2K Seed Grant 2026 &nbsp;·&nbsp; Author: Camila Duitama &nbsp;·&nbsp; Generated: {today}</p>
</header>
<main>
<div class="criteria">
  <b>Full-genome criteria applied throughout this report:</b><br>
  <code>length ≥ {MIN_LEN:,} bp</code> &nbsp;<b>AND</b>&nbsp; <code>genome_coverage ≥ {MIN_COV}</code>
  &nbsp;<b>AND</b>&nbsp; <code>is_lab_host ≠ True</code> ({n_lab} lab-adapted sequences removed)
  &nbsp;<b>AND</b>&nbsp; host is <code>Human</code> or <code>Mosquito</code>.<br>
  Other host categories (primates, rodents) and sequences with no host annotation are <b>excluded</b>.
  Data source: <b>Nextstrain DENV</b> (nextstrain.org/dengue, updated 2026-07-02), which curates
  sequences from GenBank and adds Nextclade clade/lineage annotations and per-gene coverage metrics.<br>
  Total sequences in this report: <b>{n_total:,}</b>
  ({n_human:,} human, {n_mosquito:,} mosquito).
</div>
<div class="toc"><h3>Contents</h3><ul>{toc_html}</ul></div>
{body}
</main>
<footer>Data: Nextstrain dengue (nextstrain.org) · GenBank · Generated {today}</footer>
</body>
</html>"""

out = DOCS_DIR / "index.html"
out.write_text(html, encoding="utf-8")
print(f"Site written → {out}  ({out.stat().st_size//1024} KB)")
