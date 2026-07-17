#!/usr/bin/env python3
"""
05_build_entropy_page.py
Build docs/entropy.html from data/processed/entropy/entropy_per_site.tsv.
Run after 04_entropy.py has completed.
"""

from datetime import date
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import plotly.io as pio

ENTROPY_TSV = Path("data/processed/entropy/entropy_per_site.tsv")
DOCS_DIR    = Path("docs")
DOCS_DIR.mkdir(exist_ok=True)

SERO_COLORS = {"DENV1":"#E45756","DENV2":"#F58518","DENV3":"#4C78A8","DENV4":"#72B7B2"}
HOST_COLORS = {"Human":"#4C78A8","Mosquito":"#F58518"}

GENE_ORDER = ["C","prM","E","NS1","NS2A","NS2B","NS3","NS4A","2K","NS4B","NS5"]
TARGET_GENES = {"NS4A","2K","NS4B"}

if not ENTROPY_TSV.exists():
    print(f"ERROR: {ENTROPY_TSV} not found. Run 04_entropy.py first.")
    raise SystemExit(1)

df = pd.read_csv(ENTROPY_TSV, sep="\t")
df["entropy"] = pd.to_numeric(df["entropy"], errors="coerce")

today = date.today().strftime("%d %b %Y")


def fig_html(fig, div_id):
    return pio.to_html(fig, full_html=False, include_plotlyjs=False, div_id=div_id)


def card(anchor, kind, num, title, caption, content):
    return f"""<div class="card" id="{anchor}">
  <div class="fig-label">{kind} {num} — {title}</div>
  <p class="caption">{caption}</p>
  {content}
</div>"""


# ── Figure 1: mean entropy per gene × serotype × host ─────────────────────────
def fig_summary():
    agg = (df.groupby(["serotype","gene","host"])["entropy"]
             .mean().reset_index(name="mean_H"))
    agg["gene"] = pd.Categorical(agg["gene"], categories=GENE_ORDER, ordered=True)
    agg = agg.sort_values("gene")

    fig = px.bar(agg, x="gene", y="mean_H", color="host",
                 facet_col="serotype", barmode="group",
                 color_discrete_map=HOST_COLORS,
                 labels={"mean_H":"Mean Shannon entropy (bits)",
                         "gene":"Gene","host":"Host","serotype":"Serotype"},
                 category_orders={"gene":GENE_ORDER})

    # Highlight target region
    for gene in TARGET_GENES:
        if gene in GENE_ORDER:
            fig.add_vrect(x0=GENE_ORDER.index(gene)-0.5,
                          x1=GENE_ORDER.index(gene)+0.5,
                          fillcolor="yellow", opacity=0.15, line_width=0)

    fig.update_layout(height=420, legend_title="Host", margin=dict(t=50,b=40))
    fig.for_each_annotation(lambda a: a.update(text=a.text.split("=")[-1]))
    return fig


# ── Figure 2: per-site entropy along the proteome (one panel per serotype) ────
def fig_per_site_full():
    """
    Heatmap: rows = serotype+host, columns = 20-site genomic windows,
    colour = mean entropy. Much clearer than per-site lines for showing
    which regions of the proteome are variable.
    """
    WINDOW = 20
    gene_offsets: dict[str, int] = {}
    t_off = 0
    for g in GENE_ORDER:
        gene_offsets[g] = t_off
        sub_g = df[df["gene"] == g]
        if not sub_g.empty:
            t_off += int(sub_g["site"].max())

    df2 = df.copy()
    df2["genome_site"] = df2.apply(
        lambda r: gene_offsets.get(r["gene"], 0) + r["site"], axis=1)
    df2["window"] = (df2["genome_site"] // WINDOW) * WINDOW

    # Use MM-corrected entropy if available
    ecol = "entropy_mm" if "entropy_mm" in df2.columns else "entropy"
    agg = (df2.groupby(["serotype","host","window"])[ecol]
             .mean().reset_index(name="mean_H"))
    agg["row"] = agg["serotype"] + " / " + agg["host"]

    rows_order = []
    for s in ["DENV1","DENV2","DENV3","DENV4"]:
        for h in ["Human","Mosquito"]:
            label = f"{s} / {h}"
            if label in agg["row"].values:
                rows_order.append(label)

    pivot = agg.pivot(index="row", columns="window", values="mean_H")
    pivot = pivot.reindex(rows_order)

    fig = go.Figure(go.Heatmap(
        z=pivot.values,
        x=pivot.columns.tolist(),
        y=pivot.index.tolist(),
        colorscale="RdBu_r",
        zmid=pivot.values[~np.isnan(pivot.values)].mean() if pivot.size else 1.0,
        colorbar=dict(title="Mean entropy<br>(bits)"),
        hoverongaps=False,
    ))

    # Gene labels on x-axis
    for g in GENE_ORDER:
        goff = gene_offsets.get(g, 0)
        sub_g = df[df["gene"]==g]
        if sub_g.empty: continue
        mid = goff + int(sub_g["site"].max()) / 2
        fig.add_vline(x=goff, line_dash="dot", line_color="white", line_width=1)
        if g in TARGET_GENES:
            fig.add_vrect(x0=goff, x1=goff+int(sub_g["site"].max()),
                          fillcolor="yellow", opacity=0.1, line_width=0)
        fig.add_annotation(x=mid, y=1.01, yref="paper", text=g,
                           showarrow=False, font=dict(size=9),
                           xanchor="center", yanchor="bottom")

    fig.update_layout(height=300,
                      xaxis_title=f"Genomic position ({WINDOW}-site windows)",
                      yaxis_title="",
                      margin=dict(t=40, b=40, l=140, r=20))
    return fig


# ── Figure 3: target region zoom (NS4A-2K-NS4B) ───────────────────────────────
def fig_target_zoom():
    """
    Bar chart: mean entropy ± SD per gene in the target region (NS4A-2K-NS4B),
    split by host. N is shown in the hover. Stars indicate high uncertainty (N<10).
    Uses Miller-Madow corrected entropy when available.
    """
    target_order = [g for g in GENE_ORDER if g in TARGET_GENES]
    target = df[df["gene"].isin(TARGET_GENES)].copy()
    if target.empty:
        return None

    ecol = "entropy_mm" if "entropy_mm" in target.columns else "entropy"
    agg = (target.groupby(["serotype","gene","host"])
           .agg(mean_H=(ecol,"mean"),
                sd_H=(ecol,"std"),
                n=("n_informative","median"))
           .reset_index())
    agg["gene"] = pd.Categorical(agg["gene"], categories=target_order, ordered=True)
    agg["label"] = agg["serotype"] + " — N≈" + agg["n"].round(0).astype(int).astype(str)
    agg["low_n"] = agg["n"] < 10

    serotypes = [s for s in ["DENV1","DENV2","DENV3","DENV4"]
                 if s in agg["serotype"].values]
    fig = make_subplots(rows=1, cols=len(serotypes),
                        subplot_titles=serotypes,
                        shared_yaxes=True, horizontal_spacing=0.04)

    for col, sero in enumerate(serotypes, 1):
        sub = agg[agg["serotype"]==sero].sort_values("gene")
        for host in ["Human","Mosquito"]:
            h_sub = sub[sub["host"]==host]
            if h_sub.empty: continue
            n_vals = h_sub["n"].round(0).astype(int).tolist()
            labels = [f"N≈{n}" + (" ⚠" if n<10 else "") for n in n_vals]
            fig.add_trace(go.Bar(
                x=h_sub["gene"].astype(str).tolist(),
                y=h_sub["mean_H"].tolist(),
                error_y=dict(type="data", array=h_sub["sd_H"].tolist(),
                             visible=True),
                name=host,
                legendgroup=host,
                showlegend=(col==1),
                marker_color=HOST_COLORS[host],
                text=labels,
                textposition="outside",
                textfont=dict(size=8),
                hovertemplate=(f"<b>%{{x}}</b><br>{host}<br>"
                               "Mean H = %{y:.3f} bits<br>"
                               "%{text}<extra></extra>"),
            ), row=1, col=col)

    fig.update_yaxes(title_text="Mean entropy (bits, MM-corrected)", col=1)
    fig.update_layout(height=420, barmode="group", legend_title="Host",
                      margin=dict(t=50, b=40, l=60, r=10))
    return fig


# ── Figure 4: human vs mosquito entropy scatter ───────────────────────────────
def fig_hm_scatter():
    """
    Mean entropy per gene: human vs mosquito. One point per serotype per gene.
    Stars = target genes (NS4A, 2K, NS4B). Gene labels shown next to points.
    """
    gene_means = (df.groupby(["serotype","gene","host"])["entropy"]
                  .mean().unstack("host").reset_index())
    if "Human" not in gene_means.columns or "Mosquito" not in gene_means.columns:
        return None
    gene_means = gene_means.dropna(subset=["Human","Mosquito"])
    gene_means["target"] = gene_means["gene"].isin(TARGET_GENES)
    fig = px.scatter(gene_means,
                     x="Human", y="Mosquito",
                     color="serotype", symbol="target",
                     text="gene",
                     color_discrete_map=SERO_COLORS,
                     symbol_map={True:"star", False:"circle"},
                     labels={"Human":"Mean entropy — humans (bits)",
                             "Mosquito":"Mean entropy — mosquitoes (bits)",
                             "serotype":"Serotype",
                             "target":"Target gene"},
                     hover_data={"gene":True,"Human":":.3f","Mosquito":":.3f",
                                 "serotype":True})
    fig.update_traces(textposition="top center", textfont_size=9, marker_size=10)
    vmax = max(gene_means["Human"].max(), gene_means["Mosquito"].max()) * 1.05
    fig.add_shape(type="line", x0=0, y0=0, x1=vmax, y1=vmax,
                  line=dict(dash="dash", color="grey", width=1))
    fig.update_layout(height=520, legend_title="Serotype / Target",
                      margin=dict(t=20,b=40))
    return fig


# ── Summary table ─────────────────────────────────────────────────────────────
def summary_table():
    has_std = "entropy_std" in df.columns
    agg_kw = dict(n_sites=("site","count"), n_seqs=("n_informative","median"),
                  mean_H=("entropy","mean"), max_H=("entropy","max"))
    if has_std:
        agg_kw["mean_H_std"] = ("entropy_std","mean")
    agg = (df.groupby(["serotype","gene","host"]).agg(**agg_kw)
             .reset_index().round(4))
    agg["gene"] = pd.Categorical(agg["gene"], categories=GENE_ORDER, ordered=True)
    agg = agg.sort_values(["serotype","gene","host"])
    std_th = "<th>Mean H<sub>norm</sub></th>" if has_std else ""
    tbl = ("<table><thead><tr><th>Serotype</th><th>Gene</th><th>Host</th>"
           "<th>Sites</th><th>Median seqs/site</th>"
           "<th>Mean H (bits)</th><th>Max H (bits)</th>"
           + std_th + "</tr></thead><tbody>")
    for _, r in agg.iterrows():
        hl = ' style="background:#fffde7"' if r["gene"] in TARGET_GENES else ""
        std_td = f"<td>{r['mean_H_std']:.4f}</td>" if has_std else ""
        tbl += (f"<tr{hl}><td>{r['serotype']}</td><td><b>{r['gene']}</b></td>"
                f"<td>{r['host']}</td><td>{int(r['n_sites'])}</td>"
                f"<td>{r['n_seqs']:.0f}</td>"
                f"<td>{r['mean_H']:.4f}</td><td>{r['max_H']:.4f}</td>{std_td}</tr>")
    tbl += "</tbody></table>"
    return tbl


# ── Build page ────────────────────────────────────────────────────────────────
F1 = fig_summary()
F2 = fig_per_site_full()
F3 = fig_target_zoom()
F4 = fig_hm_scatter()

body = ""
body += card("table1","Table",1,"Shannon entropy summary by serotype, gene, and host",
    "Mean entropy and number of informative sites per gene. "
    "Highlighted rows (yellow) = NS4A-2K-NS4B target region.",
    summary_table())

body += card("figure1","Figure",1,"Mean per-gene Shannon entropy by serotype and host",
    "Mean entropy across all sites per gene. Blue = human, orange = mosquito. "
    "Yellow bands mark the likely study target region (NS4A-2K-NS4B). "
    "Higher entropy = more variability at the amino acid level.",
    fig_html(F1,"f1"))

body += card("figure2","Figure",2,"Entropy heatmap across the full proteome",
    "Mean entropy per 20-site window across the polyprotein (all serotype × host combinations). "
    "Redder = more variable; bluer = more conserved. "
    "Yellow bands highlight the NS4A–2K–NS4B target region. "
    "Gene boundaries are marked with white dotted lines. "
    "Miller-Madow bias-corrected entropy is used where available.",
    fig_html(F2,"f2"))

if F3:
    body += card("figure3","Figure",3,"Mean entropy in the target region (NS4A–2K–NS4B) per serotype",
        "Bar chart of mean ± SD entropy per gene within the target region, split by host. "
        "Miller-Madow corrected entropy is used. "
        "N = approximate median number of informative sequences per site. "
        "Groups with N &lt; 10 are flagged with ⚠ — their entropy estimates have higher uncertainty "
        "(small-sample bias is partly corrected by Miller-Madow, but wide error bars indicate "
        "the mosquito data is limited).",
        fig_html(F3,"f3"))

if F4:
    body += card("figure4","Figure",4,"Human vs mosquito entropy per site",
        "Each point is one amino acid site in one gene. Points above the diagonal have higher "
        "entropy in mosquitoes; points below have higher entropy in humans. "
        "Target gene sites are shown with a different symbol.",
        fig_html(F4,"f4"))

toc = [("methods","Methods — how entropy was estimated"),
       ("table1","Table 1 — Summary statistics"),
       ("figure1","Figure 1 — Mean entropy per gene"),
       ("figure2","Figure 2 — Full proteome entropy"),
       ("figure3","Figure 3 — NS4A-2K-NS4B zoom"),
       ("figure4","Figure 4 — Human vs mosquito scatter")]
toc_html = "\n".join(f'<li><a href="#{a}">{l}</a></li>' for a,l in toc)

# Count for header
n_sero  = df["serotype"].nunique()
n_genes = df["gene"].nunique()
n_sites = len(df[df["host"]=="Human"])  # proxy

html = f"""<!DOCTYPE html><html lang="en"><head>
<meta charset="UTF-8"><meta name="viewport" content="width=device-width,initial-scale=1.0">
<title>Shannon Entropy — DENV UMR2K 2026</title>
<script src="https://cdn.plot.ly/plotly-2.35.2.min.js"></script>
<style>
:root{{--blue:#4C78A8;--bg:#f7f8fa;--card:#fff;--border:#e0e4ea;--text:#1a1a2e;--muted:#6c757d}}
*{{box-sizing:border-box;margin:0;padding:0}}
body{{font-family:'Segoe UI',Arial,sans-serif;background:var(--bg);color:var(--text);font-size:14px}}
header{{background:var(--blue);color:#fff;padding:28px 40px}}
header h1{{font-size:1.5rem;font-weight:600;margin-bottom:6px}}
header p{{opacity:.85;font-size:.88rem}}
header nav{{margin-top:12px}}
header nav a{{color:#cde4f7;text-decoration:none;margin-right:16px;font-size:.88rem}}
header nav a:hover{{text-decoration:underline}}
main{{max-width:1100px;margin:0 auto;padding:24px 20px}}
.intro{{background:#e8f4fd;border-left:4px solid var(--blue);
  padding:14px 20px;border-radius:4px;margin-bottom:20px;font-size:.9rem;line-height:1.7}}
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
table{{width:100%;border-collapse:collapse;font-size:.84rem;margin-top:8px}}
th{{background:var(--blue);color:#fff;padding:8px 12px;text-align:left;font-weight:600}}
td{{padding:6px 12px;border-bottom:1px solid var(--border)}}
tr:hover{{background:#f0f4fa}}
.methods{{background:#f0f4fa;border-left:4px solid var(--blue);
  padding:18px 24px;border-radius:4px;margin-bottom:20px;font-size:.88rem;line-height:1.7}}
.methods h2{{font-size:1rem;font-weight:600;margin-bottom:10px;color:var(--blue)}}
.methods ol{{margin-left:1.4em}}
.methods li{{margin-bottom:.6em}}
.methods ul{{margin:.4em 0 .4em 1.2em}}
.methods a{{color:var(--blue)}}
</style></head><body>
<header>
  <h1>Shannon Entropy Estimation — DENV Proteome</h1>
  <p>CNRS UMR2K Seed Grant 2026 &nbsp;·&nbsp; Author: Camila Duitama &nbsp;·&nbsp; Generated: {today}</p>
  <nav>
    <a href="index.html">← Dataset overview</a>
  </nav>
</header>
<main>

<div class="methods" id="methods">
  <h2>How entropy was estimated</h2>
  <ol>
    <li><b>Input sequences.</b>
      Full-genome DENV sequences from <b>Nextstrain</b> (nextstrain.org/dengue),
      restricted to human and mosquito hosts, ≥10,000 bp, genome coverage ≥0.95.
      Lab-grown (<code>is_lab_host = True</code>) and vaccine/chimeric strains excluded.
    </li>
    <li><b>Redundancy reduction.</b>
      Sequences were clustered at <b>98% nucleotide identity</b> using
      <a href="https://sites.google.com/view/cd-hit" target="_blank">cd-hit-est v4.8.1</a>
      (word length 8, greedy incremental algorithm).
      This reduces over-representation of densely sampled outbreaks without
      discarding biological diversity.
      Performed separately for each serotype.
    </li>
    <li><b>Translation.</b>
      Each sequence was trimmed to a multiple of three codons starting at
      nucleotide position 97 (the first codon of the C capsid protein in
      the DENV2 reference NC_001474.2 / polyprotein start).
      Translation to amino acids used
      <a href="https://emboss.sourceforge.net/" target="_blank">EMBOSS transeq v6.6.0</a>
      (<code>-frame 1 -clean</code>; ambiguous codons translated to X).
    </li>
    <li><b>Multiple sequence alignment.</b>
      All translated polyprotein sequences for a given serotype were aligned
      together using
      <a href="https://mafft.cbrc.jp/" target="_blank">MAFFT v7.526</a>
      (<code>--amino --auto</code>) with the corresponding NCBI reference
      sequence prepended as an alignment anchor
      (DENV1: NC_001477.1 | DENV2: NC_001474.2 | DENV3: NC_001475.2 | DENV4: NC_002640.1).
      This approach is equivalent to the RevTrans step used in
      Testa et al. 2026 (<i>Nat Ecol Evol</i>) for entropy purposes — both produce
      an amino-acid-level alignment; RevTrans would additionally back-project
      to a codon-aware nucleotide alignment, which is only needed for
      phylogenetic codon models (WP2, phydms).
    </li>
    <li><b>Gene boundary extraction.</b>
      The position of the reference sequence in the alignment was used to map
      nucleotide gene coordinates (from NC_001474.2) to amino acid positions
      within the polyprotein, and then to alignment columns.
      Each gene sub-alignment was extracted from the full polyprotein alignment
      by selecting the corresponding columns.
    </li>
    <li><b>Shannon entropy per site.</b>
      For each alignment column (amino acid site) in a gene:
      <ul style="margin:.5em 0 .5em 1.5em">
        <li>Gap characters (<code>-</code>), ambiguous residues (<code>X</code>),
            and stop codons (<code>*</code>) were excluded.</li>
        <li>Amino acid frequencies <i>p&#x1D43;</i> were computed from the
            remaining residues.</li>
        <li>Shannon entropy was computed as:
            <b>H = &minus;&sum; p&#x1D43; log&#x2082; p&#x1D43;</b>
            (units: bits; range 0–log&#x2082;20 &asymp; 4.32 bits).</li>
        <li>Standardised entropy was computed as:
            <b>H&#x209b;&#x209c;&#x1D00; = H / log&#x2082;(20)</b>
            (range 0–1; 0 = fully conserved, 1 = maximally diverse;
            as reported in Testa et al. 2026).</li>
        <li>Sites with fewer than 2 informative residues were excluded
            (reported as NaN).</li>
      </ul>
    </li>
    <li><b>Host comparison.</b>
      All steps above were performed independently for
      <b>human-derived</b> and <b>mosquito-derived</b> sequences,
      allowing direct comparison of site-level variability between the two host environments.
      Serotype-host groups with fewer than 5 sequences after clustering were excluded.
    </li>
  </ol>
  <p style="margin-top:10px;font-size:.82rem;color:#6c757d">
    Reference: Testa et al. (2026)
    <i>Comparative analysis of deep mutational scanning datasets in
    enteroviruses A and B identifies functional divergence and therapeutic targets.</i>
    <a href="https://doi.org/10.1038/s41559-026-02993-8" target="_blank">Nat Ecol Evol</a>.
    Code: <a href="https://github.com/QVEU/EV_DMS_Comparison" target="_blank">github.com/QVEU/EV_DMS_Comparison</a>.
  </p>
</div>
<div class="toc"><h3>Contents</h3><ul>{toc_html}</ul></div>
{body}
</main>
<footer>DENV entropy analysis · CNRS UMR2K 2026 · Generated {today}</footer>
</body></html>"""

out = DOCS_DIR / "entropy.html"
out.write_text(html, encoding="utf-8")
print(f"Entropy page → {out}  ({out.stat().st_size//1024} KB)")
