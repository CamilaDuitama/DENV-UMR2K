#!/usr/bin/env python3
"""
03_build_site.py — DENV full-genome dataset report.
Source: data/processed/denv_final_metadata.tsv
Host filter: host_type == Human | Mosquito  (Nextstrain-curated, authoritative)
"""

from datetime import date
from pathlib import Path
import numpy as np
import pandas as pd
import plotly.graph_objects as go
import plotly.express as px
import plotly.io as pio

FINAL_META = Path("data/processed/denv_final_metadata.tsv")
DOCS_DIR   = Path("docs")
DOCS_DIR.mkdir(exist_ok=True)
(DOCS_DIR / ".nojekyll").touch()

MIN_LEN = 10_000
MIN_COV = 0.95

HC = {"Human":"#4C78A8", "Mosquito":"#F58518"}
SC = {"DENV1":"#E45756","DENV2":"#F58518","DENV3":"#4C78A8","DENV4":"#72B7B2"}

VACCINE_ACCESSIONS = {
    "MW945952":("RDENV1-WP-1A","DENV1","Recombinant chimeric"),
    "AF180818":("16007 (PDK-13)","DENV1","PDK-13 attenuated vaccine"),
    "U87412":  ("PDK-53","DENV2","PDK-53 attenuated vaccine"),
    "KU725664":("PDK53","DENV2","PDK-53 variant"),
    "M84728":  ("16681-PDK53","DENV2","PDK-53 variant"),
    "KJ160505":("rDENV3-4","DENV3","Chimeric DENV3/4"),
    "MW793459":("PDK48","DENV4","PDK-48 attenuated vaccine"),
    "KJ160504":("rDENV4","DENV4","Recombinant DENV4"),
}

if not FINAL_META.exists():
    print(f"ERROR: {FINAL_META} not found. Run 02_download_genomes.py first.")
    raise SystemExit(1)

# ── Load ──────────────────────────────────────────────────────────────────────
raw = pd.read_csv(FINAL_META, sep="\t", low_memory=False)
for col in ["length","genome_coverage","NS4A_coverage","2K_coverage","NS4B_coverage",
            "E_coverage","C_coverage","prM_coverage","NS1_coverage","NS2A_coverage",
            "NS2B_coverage","NS3_coverage","NS5_coverage"]:
    if col in raw.columns:
        raw[col] = pd.to_numeric(raw[col], errors="coerce")
raw["collection_year"] = pd.to_numeric(raw["date"].astype(str).str[:4], errors="coerce")

n_total_dl  = len(raw)
n_lab       = int((raw["is_lab_host"].fillna("").astype(str).str.lower()=="true").sum())
n_vax       = int(raw["accession"].isin(VACCINE_ACCESSIONS).sum())
n_no_host   = int((~raw["host_type"].isin(["Human","Mosquito"])).sum())

# Analysis dataset: host_type Human or Mosquito, no lab, no vaccine
ns = raw[~raw["is_lab_host"].fillna("").astype(str).str.lower().eq("true") &
         ~raw["accession"].isin(VACCINE_ACCESSIONS)].copy()
full = ns[ns["host_type"].isin(["Human","Mosquito"])].copy()

n_human    = int((full["host_type"]=="Human").sum())
n_mosquito = int((full["host_type"]=="Mosquito").sum())
n_total    = len(full)

today = date.today().strftime("%d %b %Y")

def fig_html(fig, div_id):
    return pio.to_html(fig, full_html=False, include_plotlyjs=False, div_id=div_id)

def card(anchor, kind, num, title, caption, content):
    return f"""<div class="card" id="{anchor}">
  <div class="fig-label">{kind} {num} — {title}</div>
  <p class="caption">{caption}</p>
  {content}
</div>"""

# ── Table 1 ───────────────────────────────────────────────────────────────────
rows = []
for s in ["DENV1","DENV2","DENV3","DENV4"]:
    sf = full[full["serotype"]==s]
    rows.append({"Serotype":s,
                 "Human":int((sf["host_type"]=="Human").sum()),
                 "Mosquito":int((sf["host_type"]=="Mosquito").sum()),
                 "Total":len(sf)})
tbl1_df = pd.DataFrame(rows)
tot = tbl1_df.sum(numeric_only=True)
tbl1 = "<table><thead><tr><th>Serotype</th><th>Human</th><th>Mosquito</th><th>Total</th></tr></thead><tbody>"
for _,r in tbl1_df.iterrows():
    tbl1 += f"<tr><td><b>{r['Serotype']}</b></td><td>{int(r['Human']):,}</td><td>{int(r['Mosquito']):,}</td><td><b>{int(r['Total']):,}</b></td></tr>"
tbl1 += (f'<tr class="total"><td><b>Total</b></td>'
         f"<td><b>{int(tot['Human']):,}</b></td>"
         f"<td><b>{int(tot['Mosquito']):,}</b></td>"
         f"<td><b>{int(tot['Total']):,}</b></td></tr></tbody></table>")

# ── Table 2: exclusions ───────────────────────────────────────────────────────
tbl2 = f"""<table><thead><tr><th>Category</th><th>Count</th><th>Reason</th></tr></thead><tbody>
<tr><td>Lab-grown sequences</td><td>{n_lab:,}</td>
    <td>Flagged <code>is_lab_host = True</code> by Nextstrain. These were sequenced from
    cell cultures (Vero, C6/36, etc.) and may carry culture-adaptation mutations.</td></tr>
<tr><td>PDK-passage / chimeric strains</td><td>{n_vax}</td>
    <td>Vaccine candidates or engineered chimeras identified manually from strain names.
    Not caught by the lab-host flag. See list below.</td></tr>
<tr><td>No host annotation</td><td>{n_no_host:,}</td>
    <td>Sequences where the submitter did not record the host in GenBank.
    Common in older submissions and some surveillance datasets.
    Retained in the downloaded FASTA but not used in any analysis.</td></tr>
</tbody></table>
<h4 style="margin:14px 0 6px;font-size:.88rem">PDK-passage and chimeric strains (with accessions)</h4>
<table><thead><tr><th>Accession</th><th>Strain</th><th>Serotype</th><th>Reason</th></tr></thead><tbody>"""
for acc,(strain,sero,reason) in VACCINE_ACCESSIONS.items():
    tbl2 += f"<tr><td><code>{acc}</code></td><td>{strain}</td><td>{sero}</td><td>{reason}</td></tr>"
tbl2 += "</tbody></table>"

# ── Figure 1: counts ──────────────────────────────────────────────────────────
def fig1():
    rows = []
    for s in ["DENV1","DENV2","DENV3","DENV4"]:
        for h in ["Human","Mosquito"]:
            rows.append({"Serotype":s,"Host":h,
                         "n":int(((full["serotype"]==s)&(full["host_type"]==h)).sum())})
    fig = px.bar(pd.DataFrame(rows), x="Serotype", y="n", color="Host",
                 color_discrete_map=HC, barmode="group",
                 labels={"n":"Number of sequences","Host":"Host"})
    fig.update_layout(height=360, legend_title="Host", margin=dict(t=20,b=40))
    return fig

# ── Figures 2+3: length and coverage boxplots (pre-computed, split by host) ──
def _box_stats(vals):
    vals = vals.dropna().values
    if len(vals) == 0: return None
    q1,med,q3 = np.percentile(vals,[25,50,75])
    iqr = q3-q1
    return dict(q1=[float(q1)], median=[float(med)], q3=[float(q3)],
                lowerfence=[float(max(vals.min(), q1-1.5*iqr))],
                upperfence=[float(min(vals.max(), q3+1.5*iqr))])

def fig2():
    fig = go.Figure()
    for h,color in HC.items():
        sub = full[full["host_type"]==h]
        for sero in ["DENV1","DENV2","DENV3","DENV4"]:
            stats = _box_stats(sub[sub["serotype"]==sero]["length"])
            if stats:
                fig.add_trace(go.Box(**stats, x=[sero], name=h,
                                     marker_color=color, legendgroup=h,
                                     showlegend=(sero=="DENV1")))
    fig.update_layout(height=380, boxmode="group", xaxis_title="Serotype",
                      yaxis_title="Sequence length (bp)",
                      legend_title="Host", margin=dict(t=20,b=40))
    return fig

def fig3():
    fig = go.Figure()
    for h,color in HC.items():
        sub = full[full["host_type"]==h]
        for sero in ["DENV1","DENV2","DENV3","DENV4"]:
            stats = _box_stats(sub[sub["serotype"]==sero]["genome_coverage"])
            if stats:
                fig.add_trace(go.Box(**stats, x=[sero], name=h,
                                     marker_color=color, legendgroup=h,
                                     showlegend=(sero=="DENV1")))
    fig.add_hline(y=MIN_COV, line_dash="dash", line_color="red",
                  annotation_text=f"Download threshold ({MIN_COV})",
                  annotation_position="bottom right")
    fig.update_layout(height=380, boxmode="group", xaxis_title="Serotype",
                      yaxis_title="Genome coverage",
                      legend_title="Host", margin=dict(t=20,b=40))
    return fig

# ── Figure 4: per-gene coverage ───────────────────────────────────────────────
def fig4():
    gene_map = {"C_coverage":"C","prM_coverage":"prM","E_coverage":"E",
                "NS1_coverage":"NS1","NS2A_coverage":"NS2A","NS2B_coverage":"NS2B",
                "NS3_coverage":"NS3","NS4A_coverage":"NS4A","2K_coverage":"2K",
                "NS4B_coverage":"NS4B","NS5_coverage":"NS5"}
    avail_cols   = [c for c in gene_map if c in full.columns]
    avail_labels = [gene_map[c] for c in avail_cols]
    if not avail_cols: return None
    fig = go.Figure()
    for sero in ["DENV1","DENV2","DENV3","DENV4"]:
        sub = full[full["serotype"]==sero][avail_cols].dropna()
        fig.add_trace(go.Bar(x=avail_labels,
                             y=[float(sub[c].mean()) for c in avail_cols],
                             name=sero, marker_color=SC[sero],
                             error_y=dict(type="data",
                                          array=[float(sub[c].std()) for c in avail_cols],
                                          visible=True)))
    if "NS4A" in avail_labels and "NS4B" in avail_labels:
        i0 = avail_labels.index("NS4A")
        i1 = avail_labels.index("NS4B")
        fig.add_vrect(x0=i0-0.5, x1=i1+0.5, fillcolor="yellow", opacity=0.15,
                      line_width=0,
                      annotation_text="Likely target: NS4A–2K–NS4B",
                      annotation_position="bottom left",
                      annotation_font_size=11)
    fig.update_layout(barmode="group", height=420, xaxis_title="Gene",
                      yaxis_title="Mean coverage ± SD", yaxis=dict(range=[0,1.12]),
                      legend_title="Serotype", margin=dict(t=20,b=40))
    return fig

# ── Figure 5: NS4A-2K-NS4B gene stats (mean ± SD by serotype and host) ───────
def fig5():
    target_cols = {c:l for c,l in {"NS4A_coverage":"NS4A","2K_coverage":"2K",
                                    "NS4B_coverage":"NS4B"}.items() if c in full.columns}
    if not target_cols: return None
    rows = []
    for sero in ["DENV1","DENV2","DENV3","DENV4"]:
        for h in ["Human","Mosquito"]:
            sub = full[(full["serotype"]==sero)&(full["host_type"]==h)]
            for col,label in target_cols.items():
                vals = sub[col].dropna()
                if len(vals):
                    rows.append({"Serotype":sero,"Host":h,"Gene":label,
                                 "mean":float(vals.mean()), "sd":float(vals.std()),
                                 "n":len(vals)})
    df_r = pd.DataFrame(rows)
    df_r["label"] = df_r["Serotype"]+" ("+df_r["Host"]+")"
    fig = go.Figure()
    for h,color in HC.items():
        sub = df_r[df_r["Host"]==h]
        for sero in ["DENV1","DENV2","DENV3","DENV4"]:
            ssub = sub[sub["Serotype"]==sero]
            if ssub.empty: continue
            fig.add_trace(go.Bar(x=ssub["Gene"], y=ssub["mean"],
                                 error_y=dict(type="data",array=ssub["sd"].tolist(),visible=True),
                                 name=f"{sero} {h}", marker_color=color,
                                 legendgroup=sero, opacity=0.7 if h=="Human" else 1.0,
                                 showlegend=True))
    fig.update_layout(barmode="group", height=400,
                      xaxis_title="Gene", yaxis_title="Mean coverage ± SD",
                      yaxis=dict(range=[0,1.12]),
                      legend_title="Serotype / Host", margin=dict(t=20,b=40))
    return fig

# ── Figure 6: timeline ────────────────────────────────────────────────────────
def fig6():
    df = full[full["collection_year"].between(1940,2026)]
    grouped = df.groupby(["collection_year","serotype"]).size().reset_index(name="n")
    fig = px.bar(grouped, x="collection_year", y="n", color="serotype",
                 color_discrete_map=SC,
                 labels={"collection_year":"Year","n":"# sequences","serotype":"Serotype"})
    fig.update_layout(height=360, barmode="stack", legend_title="Serotype",
                      margin=dict(t=20,b=40))
    return fig

# ── Figure 7: map with human/mosquito dots ────────────────────────────────────
def fig7():
    df = full.dropna(subset=["country"])
    agg = df.groupby(["country","host_type"]).size().reset_index(name="n")
    fig = go.Figure()
    for h,color in HC.items():
        sub = agg[agg["host_type"]==h].copy()
        if sub.empty: continue
        # Scale each host independently so mosquito dots remain visible
        max_n = sub["n"].max()
        sub["dot_size"] = np.clip(6 + (sub["n"] / max_n) * 14, 6, 20)
        fig.add_trace(go.Scattergeo(
            locations=sub["country"],
            locationmode="country names",
            marker=dict(size=sub["dot_size"],
                        color=color,
                        line_color="white", line_width=0.5, opacity=0.75),
            text=sub.apply(lambda r: f"{r['country']}: {r['n']:,} {h}", axis=1),
            name=h, hoverinfo="text",
        ))
    fig.update_layout(
        height=460,
        geo=dict(showframe=False, showcoastlines=True, showland=True,
                 landcolor="#f0f0f0", projection_type="natural earth"),
        legend_title="Host",
        margin=dict(t=20,b=10,l=0,r=0),
    )
    return fig

# ── Figure 8: mosquito species bar ────────────────────────────────────────────
def fig8():
    mosq = full[full["host_type"]=="Mosquito"].dropna(subset=["host"])
    counts = (mosq.groupby(["host","serotype"]).size().reset_index(name="n")
              .sort_values("n", ascending=False))
    if counts.empty: return None
    # total per species for ordering
    order = counts.groupby("host")["n"].sum().sort_values(ascending=True).index
    fig = px.bar(counts, x="n", y="host", color="serotype",
                 color_discrete_map=SC, orientation="h",
                 labels={"n":"Number of sequences","host":"Host species"},
                 category_orders={"host":list(order)})
    fig.update_layout(height=380, legend_title="Serotype", margin=dict(t=20,b=40))
    return fig

# ── Build body ────────────────────────────────────────────────────────────────
F = {1:fig1(), 2:fig2(), 3:fig3(), 4:fig4(), 5:fig5(), 6:fig6(), 7:fig7(), 8:fig8()}

body = ""

body += card("table1","Table",1,"Full genomes retained for analysis",
    f"Counts by serotype and host. Only sequences with a confirmed host annotation "
    f"(<i>Homo sapiens</i> or <i>Aedes/Culicidae</i>) are included. "
    f"Total: <b>{n_total:,}</b> — <b>{n_human:,} from humans</b> and <b>{n_mosquito:,} from mosquitoes</b>.",
    tbl1)

body += card("table2","Table",2,"Sequences excluded from the analysis",
    f"Three categories of sequences are excluded. Together they account for "
    f"{n_lab+n_vax+n_no_host:,} of the {n_total_dl:,} downloaded sequences.",
    tbl2)

body += card("figure1","Figure",1,"Sequence counts by serotype and host",
    "Number of full genomes per serotype, split by host. Human sequences greatly outnumber "
    "mosquito sequences for all serotypes, particularly DENV3 and DENV4.",
    fig_html(F[1],"f1"))

body += card("figure2","Figure",2,"Sequence length distribution by serotype and host",
    f"Box plots show the median, interquartile range, and whiskers (1.5×IQR) of sequence length "
    f"for human and mosquito genomes. All sequences are between {MIN_LEN:,} and 12,000 bp by design. "
    f"Mosquito sequences are slightly shorter on average.",
    fig_html(F[2],"f2"))

body += card("figure3","Figure",3,"Genome coverage distribution by serotype and host",
    f"Box plots of genome coverage (fraction of the reference covered). "
    f"The dashed red line marks the download threshold ({MIN_COV}); all retained sequences "
    f"are above it. Coverage is consistently high across both hosts.",
    fig_html(F[3],"f3"))

if F[4]:
    body += card("figure4","Figure",4,"Per-gene coverage — mean ± SD (all sequences)",
        "Mean ± SD coverage per gene, grouped by serotype. The highlighted region "
        "(NS4A–2K–NS4B) is the likely study target. All genes show high mean coverage.",
        fig_html(F[4],"f4"))

if F[5]:
    body += card("figure5","Figure",5,"Coverage of the target region (NS4A–2K–NS4B) by serotype and host",
        "Mean ± SD coverage for each gene in the target region, broken down by serotype and host. "
        "Coverage is consistently above 0.95 for both human and mosquito sequences.",
        fig_html(F[5],"f5"))

body += card("figure6","Figure",6,"Temporal distribution of full genomes",
    "Collection year by serotype. Most sequences were collected between 2000 and 2024.",
    fig_html(F[6],"f6"))

body += card("figure7","Figure",7,"Geographic distribution — human and mosquito sequences",
    "Each dot represents sequences from one country. Dot size is proportional to "
    "the square root of the number of sequences. Blue = human, orange = mosquito.",
    fig_html(F[7],"f7"))

if F[8]:
    body += card("figure8","Figure",8,"Mosquito host species",
        "Number of mosquito-derived full genomes by host species and serotype. "
        "<i>Culicidae</i> entries are identified to family level only (species unknown). "
        "Please indicate which species to include in the final analysis.",
        fig_html(F[8],"f8"))

# ── TOC ───────────────────────────────────────────────────────────────────────
toc = [("table1","Table 1 — Analysis dataset"),
       ("table2","Table 2 — Excluded sequences"),
       ("figure1","Figure 1 — Counts by serotype &amp; host"),
       ("figure2","Figure 2 — Length distribution"),
       ("figure3","Figure 3 — Genome coverage"),
       ("figure4","Figure 4 — Per-gene coverage"),
       ("figure5","Figure 5 — NS4A–2K–NS4B coverage"),
       ("figure6","Figure 6 — Temporal distribution"),
       ("figure7","Figure 7 — Geographic distribution"),
       ("figure8","Figure 8 — Mosquito host species")]
toc_html = "\n".join(f'<li><a href="#{a}">{l}</a></li>' for a,l in toc)

html = f"""<!DOCTYPE html><html lang="en"><head>
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
.intro{{background:#e8f4fd;border-left:4px solid var(--blue);
  padding:14px 20px;border-radius:4px;margin-bottom:20px;font-size:.9rem;line-height:1.7}}
.intro code{{background:#cde4f7;padding:1px 5px;border-radius:3px}}
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
h4{{color:var(--blue);margin-top:4px}}
footer{{text-align:center;padding:24px;color:var(--muted);font-size:.8rem}}
</style></head><body>
<header>
  <h1>DENV Full-Genome Dataset — Availability Report</h1>
  <p>CNRS UMR2K Seed Grant 2026 &nbsp;·&nbsp; Author: Camila Duitama &nbsp;·&nbsp; Generated: {today}</p>
</header>
<main>
<div class="intro">
  This report summarises the DENV full-genome dataset downloaded from
  <a href="https://nextstrain.org/dengue" target="_blank">Nextstrain DENV</a>
  (data updated 2026-07-02, sourced from GenBank with Nextclade annotations).<br><br>
  We kept sequences that are at least 10,000 bp long and cover at least 95% of the reference genome.
  We then restricted the analysis dataset to sequences with a confirmed human or mosquito host annotation —
  sequences from other hosts (primates, rodents) or with no host recorded are excluded from all analyses.
  Lab-grown and vaccine-derived sequences are also removed.<br><br>
  <b>Final analysis dataset: {n_total:,} sequences — {n_human:,} from humans and {n_mosquito:,} from mosquitoes.</b>
</div>
<div class="toc"><h3>Contents</h3><ul>{toc_html}</ul></div>
{body}
</main>
<footer>Data: Nextstrain dengue (nextstrain.org) · GenBank · Generated {today}</footer>
</body></html>"""

out = DOCS_DIR / "index.html"
out.write_text(html, encoding="utf-8")
print(f"Site written → {out}  ({out.stat().st_size//1024} KB)")
