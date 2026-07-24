#!/usr/bin/env python3
"""
05_build_entropy_page.py
Build docs/entropy.html from data/processed/entropy/entropy_per_site.tsv.
Run after 04_entropy.py has completed.
"""

from datetime import date
from pathlib import Path

import math
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import plotly.io as pio

PIPELINE_SVG = """<svg xmlns="http://www.w3.org/2000/svg" width="580" height="895"
     font-family="Arial, Helvetica, sans-serif">
  <defs>
    <marker id="arr" markerWidth="10" markerHeight="7" refX="9" refY="3.5" orient="auto">
      <polygon points="0 0,10 3.5,0 7" fill="#555"/>
    </marker>
    <marker id="arr-red" markerWidth="10" markerHeight="7" refX="9" refY="3.5" orient="auto">
      <polygon points="0 0,10 3.5,0 7" fill="#c0392b"/>
    </marker>
  </defs>

  <!-- INPUT  y=10..62 -->
  <rect x="60" y="10" width="460" height="52" rx="8" fill="#f5a623" stroke="#c47c00" stroke-width="2"/>
  <text x="290" y="29" text-anchor="middle" font-weight="bold" font-size="13">INPUT</text>
  <text x="290" y="48" text-anchor="middle" font-size="12">Nextstrain DENV sequences (full genomes)</text>

  <line x1="290" y1="62" x2="290" y2="83" stroke="#555" stroke-width="2" marker-end="url(#arr)"/>

  <!-- FILTER  y=83..149 -->
  <rect x="60" y="83" width="460" height="66" rx="8" fill="#4a90d9" stroke="#2c5f8a" stroke-width="2"/>
  <text x="290" y="102" text-anchor="middle" font-weight="bold" font-size="13" fill="white">FILTER</text>
  <text x="290" y="120" text-anchor="middle" font-size="12" fill="white">Human / mosquito host &#x00B7; length &#x2265;10&#x202F;kbp &#x00B7; coverage &#x2265;0.95</text>
  <text x="290" y="137" text-anchor="middle" font-size="12" fill="white">Vaccine strains and lab-grown sequences excluded</text>

  <line x1="290" y1="149" x2="290" y2="170" stroke="#555" stroke-width="2" marker-end="url(#arr)"/>

  <!-- SPLIT  y=170..222 -->
  <rect x="60" y="170" width="460" height="52" rx="8" fill="#4a90d9" stroke="#2c5f8a" stroke-width="2"/>
  <text x="290" y="189" text-anchor="middle" font-weight="bold" font-size="13" fill="white">SPLIT</text>
  <text x="290" y="207" text-anchor="middle" font-size="12" fill="white">Independently by serotype (DENV1&#x2013;4) and host (Human / Mosquito)</text>

  <line x1="290" y1="222" x2="290" y2="246" stroke="#888" stroke-width="2" stroke-dasharray="5,3" marker-end="url(#arr)"/>

  <!-- LOOP border  y=248..783 -->
  <rect x="38" y="248" width="504" height="535" rx="10" fill="none" stroke="#bbb" stroke-width="1.5" stroke-dasharray="7,4"/>
  <text x="290" y="264" text-anchor="middle" font-size="10" fill="#999" font-style="italic">repeated for each serotype &#x00D7; host combination</text>

  <!-- CLUSTER  y=270..336 -->
  <rect x="80" y="270" width="420" height="66" rx="8" fill="#4a90d9" stroke="#2c5f8a" stroke-width="2"/>
  <text x="290" y="289" text-anchor="middle" font-weight="bold" font-size="13" fill="white">CLUSTER</text>
  <text x="290" y="307" text-anchor="middle" font-size="12" fill="white">Cluster at 98% nucleotide identity (cd-hit-est)</text>
  <text x="290" y="324" text-anchor="middle" font-size="12" fill="white">Reduces outbreak over-representation &#x00B7; per serotype</text>

  <line x1="290" y1="336" x2="290" y2="360" stroke="#555" stroke-width="2" marker-end="url(#arr)"/>

  <!-- DECISION diamond  center=(290,390)  half-w=108 half-h=30 -->
  <polygon points="290,360 398,390 290,420 182,390" fill="#e74c3c" stroke="#a93226" stroke-width="2"/>
  <text x="290" y="385" text-anchor="middle" font-weight="bold" font-size="12" fill="white">n &lt; 5 seqs?</text>
  <text x="290" y="400" text-anchor="middle" font-size="10" fill="white">after clustering</text>

  <!-- YES branch -->
  <line x1="398" y1="390" x2="422" y2="390" stroke="#c0392b" stroke-width="2" marker-end="url(#arr-red)"/>
  <text x="410" y="382" text-anchor="middle" font-size="10" fill="#c0392b" font-weight="bold">Yes</text>
  <rect x="424" y="370" width="88" height="40" rx="6" fill="#95a5a6" stroke="#717d7e" stroke-width="1.5"/>
  <text x="468" y="386" text-anchor="middle" font-size="11" fill="white" font-weight="bold">EXCLUDE</text>
  <text x="468" y="401" text-anchor="middle" font-size="10" fill="white">group skipped</text>

  <!-- NO branch -->
  <line x1="290" y1="420" x2="290" y2="444" stroke="#555" stroke-width="2" marker-end="url(#arr)"/>
  <text x="302" y="437" font-size="10" fill="#333" font-weight="bold">No</text>

  <!-- TRANSLATE  y=444..510 -->
  <rect x="80" y="444" width="420" height="66" rx="8" fill="#4a90d9" stroke="#2c5f8a" stroke-width="2"/>
  <text x="290" y="463" text-anchor="middle" font-weight="bold" font-size="13" fill="white">TRANSLATE</text>
  <text x="290" y="481" text-anchor="middle" font-size="12" fill="white">Translate polyprotein to amino acids</text>
  <text x="290" y="498" text-anchor="middle" font-size="12" fill="white">Start: nucleotide 97 &#x00B7; ambiguous residues &#x2192; X</text>

  <line x1="290" y1="510" x2="290" y2="531" stroke="#555" stroke-width="2" marker-end="url(#arr)"/>

  <!-- ALIGN  y=531..597 -->
  <rect x="80" y="531" width="420" height="66" rx="8" fill="#4a90d9" stroke="#2c5f8a" stroke-width="2"/>
  <text x="290" y="550" text-anchor="middle" font-weight="bold" font-size="13" fill="white">ALIGN</text>
  <text x="290" y="568" text-anchor="middle" font-size="12" fill="white">Multiple amino acid alignment (MAFFT --amino)</text>
  <text x="290" y="585" text-anchor="middle" font-size="12" fill="white">NCBI reference sequence prepended as anchor</text>

  <line x1="290" y1="597" x2="290" y2="618" stroke="#555" stroke-width="2" marker-end="url(#arr)"/>

  <!-- EXTRACT GENES  y=618..684 -->
  <rect x="80" y="618" width="420" height="66" rx="8" fill="#4a90d9" stroke="#2c5f8a" stroke-width="2"/>
  <text x="290" y="637" text-anchor="middle" font-weight="bold" font-size="13" fill="white">EXTRACT GENES</text>
  <text x="290" y="655" text-anchor="middle" font-size="12" fill="white">Map gene coordinates to alignment columns</text>
  <text x="290" y="672" text-anchor="middle" font-size="12" fill="white">Extract one sub-alignment per gene</text>

  <line x1="290" y1="684" x2="290" y2="705" stroke="#555" stroke-width="2" marker-end="url(#arr)"/>

  <!-- ENTROPY  y=705..775 -->
  <rect x="80" y="705" width="420" height="70" rx="8" fill="#1a5276" stroke="#154360" stroke-width="2"/>
  <text x="290" y="724" text-anchor="middle" font-weight="bold" font-size="13" fill="white">ENTROPY</text>
  <text x="290" y="742" text-anchor="middle" font-size="12" fill="white">H = &#x2212;&#x2211;&#x202F;p&#x202F;log&#x2082;&#x202F;p (bits, range 0&#x2013;4.32)</text>
  <text x="290" y="759" text-anchor="middle" font-size="12" fill="white">Gaps, ambiguous (X), and stop (*) residues excluded</text>

  <!-- end loop -->
  <line x1="290" y1="775" x2="290" y2="796" stroke="#555" stroke-width="2" marker-end="url(#arr)"/>

  <!-- OUTPUT  y=796..860 -->
  <rect x="60" y="796" width="460" height="64" rx="8" fill="#27ae60" stroke="#1a6e3c" stroke-width="2"/>
  <text x="290" y="815" text-anchor="middle" font-weight="bold" font-size="13" fill="white">OUTPUT</text>
  <text x="290" y="833" text-anchor="middle" font-size="12" fill="white">Per-site entropy table (H &#x00B7; H<tspan dy="4" font-size="9">std</tspan><tspan dy="-4">)</tspan></text>
  <text x="290" y="850" text-anchor="middle" font-size="12" fill="white">Per gene &#x00B7; per serotype &#x00B7; per host</text>

  <!-- Legend  y=876 -->
  <g transform="translate(60,876)">
    <rect x="0" y="0" width="12" height="12" rx="2" fill="#f5a623"/>
    <text x="16" y="11" font-size="10" fill="#555">Input / Output</text>
    <rect x="112" y="0" width="12" height="12" rx="2" fill="#4a90d9"/>
    <text x="128" y="11" font-size="10" fill="#555">Processing step</text>
    <polygon points="232,6 246,0 260,6 246,12" fill="#e74c3c"/>
    <text x="264" y="11" font-size="10" fill="#555">Decision gate</text>
    <rect x="358" y="0" width="12" height="12" rx="2" fill="#95a5a6"/>
    <text x="374" y="11" font-size="10" fill="#555">Excluded</text>
  </g>
</svg>"""


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


def card(anchor, kind, num, title, caption, content, collapsible=False, finding=None):
    if collapsible:
        body = (f'<details><summary style="cursor:pointer;color:var(--blue);font-size:.85rem;margin-bottom:6px">Show / hide table</summary>\n'
                f'  <p class="caption">{caption}</p>\n  {content}\n</details>')
    else:
        body = f'  <p class="caption">{caption}</p>\n  {content}'
    if finding:
        body += f'\n  <div class="finding"><span class="finding-label">Key finding</span>{finding}</div>'
    return f"""<div class="card" id="{anchor}">
  <div class="fig-label">{kind} {num} — {title}</div>
{body}
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

    # Add genome-wide mean per serotype as a dashed reference line
    sero_order = [s for s in ["DENV1","DENV2","DENV3","DENV4"]
                  if s in df["serotype"].unique()]
    # Gene-unweighted mean: mean of per-gene means (appropriate for this bar chart
    # where each bar = one gene regardless of length).
    # Site-weighted mean would be dominated by long, conserved genes (e.g. NS3 in DENV2).
    gw_h = (df[df["host"]=="Human"].groupby(["serotype","gene"])["entropy"]
            .mean().groupby("serotype").mean())
    gw_m = (df[df["host"]=="Mosquito"].groupby(["serotype","gene"])["entropy"]
            .mean().groupby("serotype").mean())
    for col_i, sero in enumerate(sero_order, 1):
        if sero in gw_h.index:
            fig.add_hline(
                y=gw_h[sero], row=1, col=col_i,
                line_dash="dot", line_color=HOST_COLORS["Human"], line_width=1.5,
                annotation_text=f"human mean {gw_h[sero]:.2f}",
                annotation_font_size=7,
                annotation_position="top right",
            )
        if sero in gw_m.index:
            fig.add_hline(
                y=gw_m[sero], row=1, col=col_i,
                line_dash="dash", line_color=HOST_COLORS["Mosquito"], line_width=1.2,
                annotation_text=f"mosq. mean {gw_m[sero]:.2f}",
                annotation_font_size=7,
                annotation_position="bottom right",
            )
    fig.update_layout(height=420, legend_title="Host", margin=dict(t=50,b=60))
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


# ── Figure 3: target region dot plot (NS4A-2K-NS4B) ──────────────────────────
LOG2_20 = math.log2(20)

def fig_target_zoom():
    """
    Cleveland dot plot: mean entropy per gene in target region (NS4A-2K-NS4B),
    one dot per serotype, colored by serotype. Human sequences only, n>=10.
    Vertical line = genome-wide mean (human). Error bars = SD.
    """
    target_order = [g for g in GENE_ORDER if g in TARGET_GENES]
    target = df[(df["gene"].isin(TARGET_GENES)) & (df["host"]=="Human")].copy()
    if target.empty:
        return None

    ecol = "entropy_mm" if "entropy_mm" in target.columns else "entropy"
    agg = (target.groupby(["serotype","gene"])
           .agg(mean_H=(ecol,"mean"), sd_H=(ecol,"std"), n=("n_informative","median"))
           .reset_index())
    agg = agg[agg["n"] >= 10]
    if agg.empty:
        return None
    agg["gene"] = pd.Categorical(agg["gene"], categories=target_order, ordered=True)

    overall_gw = df[df["host"]=="Human"][ecol].mean()

    fig = go.Figure()
    serotypes = [s for s in ["DENV1","DENV2","DENV3","DENV4"]
                 if s in agg["serotype"].values]
    # y-offsets so dots for different serotypes don't overlap
    y_offsets = {s: off for s, off in
                 zip(serotypes, [-0.28, -0.09, 0.09, 0.28][:len(serotypes)])}
    for sero in serotypes:
        sub = agg[agg["serotype"]==sero].sort_values("gene")
        y_vals = [target_order.index(g) + y_offsets[sero]
                  for g in sub["gene"].astype(str)]
        fig.add_trace(go.Scatter(
            x=sub["mean_H"].tolist(),
            y=y_vals,
            mode="markers",
            marker=dict(size=11, color=SERO_COLORS[sero]),
            error_x=dict(type="data", array=sub["sd_H"].fillna(0).tolist(), visible=True,
                         color=SERO_COLORS[sero], thickness=1.5, width=4),
            name=sero,
            hovertemplate=(f"<b>%{{customdata[0]}}</b><br>{sero}<br>"
                           "Mean H = %{x:.3f} bits<extra></extra>"),
            customdata=[[g] for g in sub["gene"].astype(str)],
        ))

    fig.add_vline(x=overall_gw, line_dash="dot", line_color="grey", line_width=1.5,
                  annotation_text=f"genome-wide mean<br>(human, {overall_gw:.2f} bits)",
                  annotation_font_size=9, annotation_position="top left")
    fig.add_vline(x=0, line_color="white", line_width=0)  # force axis from 0
    fig.update_xaxes(title_text="Mean entropy H (bits)", range=[0, LOG2_20 * 0.45],
                     tickvals=[0, 0.5, 1.0, 1.5], ticktext=["0<br><i>conserved</i>","0.5","1.0","1.5"])
    fig.update_yaxes(
        tickvals=list(range(len(target_order))),
        ticktext=target_order,
        tickfont=dict(size=13),
        range=[-0.5, len(target_order) - 0.5],
    )
    fig.update_layout(height=300, legend_title="Serotype",
                      margin=dict(t=30, b=50, l=80, r=20))
    return fig


# ── Figure 4: human vs mosquito entropy scatter ───────────────────────────────
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
        n = r["n_seqs"]
        if n < 5:
            # Unreliable: too few sequences per site — strikethrough, dark grey
            style = "background:#f5f5f5;color:#aaa;text-decoration:line-through"
            flag = " &#x26A0;"
        elif n < 10:
            # Caution: low N, entropy likely underestimated
            style = "background:#fff3cd"
            flag = " &#x26A0;"
        elif r["gene"] in TARGET_GENES:
            style = "background:#fffde7"
            flag = ""
        else:
            style = ""
            flag = ""
        std_td = f"<td>{r['mean_H_std']:.4f}</td>" if has_std else ""
        tbl += (f'<tr style="{style}"><td>{r["serotype"]}</td><td><b>{r["gene"]}</b></td>'
                f'<td>{r["host"]}</td><td>{int(r["n_sites"])}</td>'
                f'<td>{r["n_seqs"]:.0f}{flag}</td>'
                f'<td>{r["mean_H"]:.4f}</td><td>{r["max_H"]:.4f}</td>{std_td}</tr>')
    tbl += ("</tbody></table>"
            '<p style="margin-top:8px;font-size:.8rem;color:#555">'
            "&#x26A0; = median informative sequences per site &lt;&thinsp;10 "
            "(entropy estimates unreliable; use with caution). "
            "Strikethrough = median &lt;&thinsp;5 (essentially meaningless). "
            "Yellow = NS4A&ndash;2K&ndash;NS4B target region.</p>")
    return tbl



# ── Build page ────────────────────────────────────────────────────────────────
F1 = fig_summary()
F2 = fig_per_site_full()
F3 = fig_target_zoom()

body = ""

body += card("figure1","Figure",1,"Mean per-gene Shannon entropy by serotype and host",
    "Mean Shannon entropy H (bits) per gene, split by serotype and host. "
    "Dotted line = gene-unweighted mean entropy for human sequences (mean of the 11 gene means, per serotype). "
    "Dashed line = genome-wide mean for mosquito sequences (per serotype). "
    "Yellow bands = NS4A–2K–NS4B target region.",
    fig_html(F1,"f1"),
    finding=("Most variable genes (human): <b>prM</b> (1.14 bits), <b>E</b> (1.10), <b>NS2A</b> (1.10), <b>NS4B</b> (1.08) — "
             "genes under direct antibody or immune selection. "
             "Most conserved: <b>2K</b> (0.93), <b>NS2B</b> (1.02), <b>NS4A</b> (1.02). "
             "Note: all genes cluster within a narrow range (0.93–1.14 bits); differences are real but modest. "
             "Target region (NS4A–2K–NS4B) mean = 1.01 bits, slightly <em>below</em> the genome-wide mean (1.06 bits)."))

body += card("figure2","Figure",2,"Entropy heatmap across the full proteome",
    "Mean H per 20-site sliding window across the DENV polyprotein (all serotype × host combinations). "
    "Redder = more variable; bluer = more conserved. Yellow bands = target region (NS4A–2K–NS4B).",
    fig_html(F2,"f2"),
    finding=("<b>DENV2</b> is notably more conserved (0.77\u202fbits) than DENV1/3/4 (1.1\u20131.2\u202fbits) \u2014 visible as the consistently blue row across all panels. "
             "The target region shows a relatively conserved band across serotypes."))

if F3:
    body += card("figure3","Figure",3,"Mean entropy in the target region (NS4A–2K–NS4B) per serotype",
        "Mean H ± SD per gene in the NS4A–2K–NS4B target region (human sequences only). "
        "Serotype groups with fewer than 10 informative sequences are excluded. "
        "Dashed reference lines = genome-wide mean H per serotype.",
        fig_html(F3,"f3"),
        finding=("Conservation within the target region <b>varies by serotype</b>: "
                 "in DENV1/3/4, <b>2K</b> (signal peptide) is the most conserved gene (lowest H). "
                 "In DENV2, NS4A and NS4B are more conserved than 2K (0.80 and 0.81 vs 1.04 bits). "
                 "<b>DENV2 is the most conserved serotype</b> across NS4A and NS4B, "
                 "but its 2K segment is an exception: 2K has <em>higher</em> entropy (1.04 bits) "
                 "than NS4A (0.80) and NS4B (0.81) in DENV2 only. "
                 "This is verified from the raw alignment (620 sequences, ~565 informative per site) "
                 "and reflects genuine variability at 18 of 21 sites. "
                 "Three sites in DENV2 2K are highly conserved (G/L/Q at >97%), "
                 "but the remaining 18 carry 2–4 amino acids each at substantial frequency. "
                 "Only human sequences shown (mosquito N < 10 after clustering)."))


body += card("table1","Table",1,"Full entropy statistics by serotype, gene, and host",
    "Mean H (bits), max H, and H<sub>norm</sub> = H\u202f/\u202flog\u2082(20) (range 0\u20131) per gene \u00d7 host \u00d7 serotype. "
    "Rows are colour-coded by data reliability: "
    "white = reliable (median seqs/site \u226510); "
    "amber = caution (5\u20139 seqs/site, entropy likely underestimated); "
    "grey strikethrough = unreliable (&lt;5 seqs/site). "
    "Yellow = NS4A\u20132K\u2013NS4B target region.",
    summary_table(), collapsible=True)

toc = [("methods","Methods — pipeline"),
       ("figure1","Figure 1 — Mean entropy per gene"),
       ("figure2","Figure 2 — Full proteome entropy"),
       ("figure3","Figure 3 — NS4A-2K-NS4B zoom"),
       ("table1","Table 1 — Full statistics")]
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
.pipeline-svg{{max-width:100%;display:block;margin:14px auto}}
.finding{{background:#e8f4fd;border-left:3px solid var(--blue);padding:8px 14px;margin-top:10px;border-radius:0 4px 4px 0;font-size:.88rem;line-height:1.55}}
.finding-label{{font-weight:700;color:var(--blue);margin-right:5px}}
</style></head><body>
<header>
  <h1>Shannon Entropy Estimation — DENV Proteome</h1>
  <p>CNRS UMR2K Seed Grant 2026 &nbsp;·&nbsp; Author: Camila Duitama &nbsp;·&nbsp; Generated: {today}</p>
  <nav>
    <a href="index.html">← Dataset overview</a>
    &nbsp;&middot;&nbsp;
    <a href="https://camiladuitama.github.io/DENV-UMR2K/entropy.html" target="_blank">&#127760; Share this page</a>
  </nav>
</header>
<main>
<div class="toc"><h3>Contents</h3><ul>{toc_html}</ul></div>

<div class="methods" id="methods">
  <h2>Pipeline</h2>
  <div class="pipeline-svg">{PIPELINE_SVG}</div>
  </div>
  <p style="margin-top:8px;font-size:.82rem;color:#6c757d">
    H in bits (0 = conserved &rarr; 4.32 = maximally diverse). H<sub>std</sub> = H&thinsp;/&thinsp;log&#x2082;20 (0&ndash;1).
    Analysed independently per serotype &times; host. Min 5 sequences per group after clustering.
    Approach adapted from <a href="https://doi.org/10.1038/s41559-026-02993-8" target="_blank">Testa et al. 2026 (<i>Nat Ecol Evol</i>)</a>.
    <b>Key difference from Testa et al.:</b> we align at the amino acid level (MAFFT <code>--amino</code>)
    and compute entropy directly, skipping the RevTrans step (back-translation to a codon alignment).
    RevTrans is only needed for phylogenetic codon models (phydms, future WP2) and is not required for amino acid entropy.
  </p>
</div>
{body}
</main>
<footer>DENV entropy analysis · CNRS UMR2K 2026 · Generated {today}</footer>
</body></html>"""

out = DOCS_DIR / "entropy.html"
out.write_text(html, encoding="utf-8")
print(f"Entropy page → {out}  ({out.stat().st_size//1024} KB)")
