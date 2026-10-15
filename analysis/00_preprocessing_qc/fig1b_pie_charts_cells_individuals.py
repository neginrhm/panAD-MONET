"""
Fig. 1b pie charts: atlas composition by major cell type and donor counts,
drawn side by side as (nuclei | individuals).

Input (working directory):
  13ad_3m_5k_atlas_final_renamed_endo_meta_filled_no_bc.h5ad
Output dir:
  pie_charts_atlas_final_renamed_with_individuals/
"""

import h5py
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib import rcParams
import os

rcParams.update({
    "font.family":     "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "pdf.fonttype":    42,
    "svg.fonttype":    "none",
})

H5AD = "13ad_3m_5k_atlas_final_renamed_endo_meta_filled_no_bc.h5ad"
OUT_DIR = "pie_charts_atlas_final_renamed_with_individuals"
os.makedirs(OUT_DIR, exist_ok=True)


# ── I/O helper ────────────────────────────────────────────────────────────────
def read_col(f, key):
    grp = f["obs"][key]
    if "categories" in grp:
        cats = [c.decode() if isinstance(c, bytes) else c for c in grp["categories"][:]]
        codes = grp["codes"][:]
        return pd.Categorical.from_codes(codes, cats)
    vals = grp[:]
    if vals.dtype.kind in ("S", "O"):
        vals = [v.decode() if isinstance(v, bytes) else v for v in vals]
    return vals


# ── load ──────────────────────────────────────────────────────────────────────
print("Reading h5ad …")
cols_needed = ["Celltype_markers", "Study", "Sex", "Tissue",
               "Tissue_standardized", "Diagnosis", "dx",
               "external_donor_name_label"]
with h5py.File(H5AD, "r") as f:
    data = {k: read_col(f, k) for k in cols_needed}

obs = pd.DataFrame(data)
obs = obs.rename(columns={"external_donor_name_label": "donor"})
print(f"Loaded {len(obs):,} cells, {obs['donor'].nunique():,} donors")


# ── colour palettes (match originals) ─────────────────────────────────────────
COLORS = {
    # Celltype_markers
    "Excitatory":      "#d62728",
    "Oligodendrocyte": "#2ca02c",
    "Inhibitory":      "#1f77b4",
    "Astrocyte":       "#ff7f0e",
    "Microglia":       "#9467bd",
    "OPC":             "#8c564b",
    "Endothelial":     "#e377c2",

    # Study
    "ROSMAP":               "#1f77b4",
    "SEA-AD":               "#2ca02c",
    "ADDS":                 "#d62728",
    "AD_Leng_2021":         "#ff7f0e",
    "AD_Zhou_2020":         "#9467bd",
    "ADvas_Yang_2021":      "#8c564b",
    "AD_Mathys_2019":       "#aec7e8",
    "AD_Swarup":            "#98df8a",
    "Morabito_2020":        "#ffbb78",
    "AD_Otero-Garcia_2020": "#ff9896",
    "AD_Grubman_2019":      "#c5b0d5",
    "Olah_microglia_AD":    "#c7e9c0",
    "AD_Lau_2020":          "#b0e0e6",

    # Sex
    "Female":  "#e377c2",
    "Male":    "#aec7e8",

    # Tissue (raw)
    "dorsolateral prefrontal cortex": "#9467bd",
    "middle temporal gyrus":          "#1f77b4",
    "pfc":                            "#2ca02c",
    "pcc":                            "#ff7f0e",
    "fcx":                            "#d62728",
    "Superior frontal gyrus":         "#c5b0d5",
    "Hippocampus":                    "#7f7f7f",
    "unknown":                        "#c7c7c7",
    "Entorhinal cortex":              "#17becf",
    "Entorhinal Cortex":              "#bcbd22",

    # Tissue_standardized
    "Frontal cortex":               "#1f77b4",
    "Temporal cortex":              "#2ca02c",
    "Posterior cingulate cortex":   "#ff7f0e",
    "Hippocampus_std":              "#8c564b",   # keyed separately below
    "Entorhinal cortex_std":        "#9467bd",

    # Diagnosis (fine)
    "AD":                    "#d62728",
    "Control":               "#2ca02c",
    "NCI":                   "#17becf",
    "MCI":                   "#ff7f0e",
    "Dementia":              "#e377c2",
    "No Dementia":           "#bcbd22",
    "DSAD":                  "#9467bd",
    "AD, Vascular Dementia": "#ffbb78",
    "PD":                    "#c7c7c7",
    "Vascular Dementia":     "#8c564b",

    # dx (broad)
    "Dementia_dx": "#d62728",
    "Control_dx":  "#2ca02c",
    "MCI_dx":      "#ff7f0e",
}

# per-column colour getter with fallback tab20
_TAB20 = plt.cm.tab20.colors

def get_colors(labels, col_key=None):
    out = []
    for i, lab in enumerate(labels):
        if lab in COLORS:
            out.append(COLORS[lab])
        elif col_key == "Tissue_standardized" and lab in ("Hippocampus", "Entorhinal cortex"):
            out.append(COLORS.get(lab + "_std", _TAB20[i % 20]))
        elif col_key == "dx":
            out.append(COLORS.get(lab + "_dx", _TAB20[i % 20]))
        else:
            out.append(_TAB20[i % 20])
    return out


# ── core pie function ──────────────────────────────────────────────────────────
def make_pie(ax, counts, title, subtitle, colors, pct_threshold=2.0):
    """Draw one pie; labels on wedge only when large enough."""
    labels = list(counts.index)
    sizes  = list(counts.values)
    total  = sum(sizes)

    def autopct(pct):
        return f"{pct:.1f}%" if pct >= pct_threshold else ""

    wedges, texts, autotexts = ax.pie(
        sizes,
        labels=None,
        colors=colors,
        autopct=autopct,
        startangle=90,
        pctdistance=0.78,
        wedgeprops=dict(linewidth=0.8, edgecolor="white"),
        textprops=dict(fontsize=10),
    )
    for at in autotexts:
        at.set_fontsize(10)

    ax.set_title(f"{title}\n({subtitle})", fontsize=13, pad=10)

    legend_labels = [f"{lab}  ({cnt:,}  {cnt/total*100:.1f}%)"
                     for lab, cnt in zip(labels, sizes)]
    patches = [mpatches.Patch(color=c, label=l)
               for c, l in zip(colors, legend_labels)]
    ax.legend(handles=patches,
              loc="center left",
              bbox_to_anchor=(1.02, 0.5),
              frameon=False,
              fontsize=9,
              title_fontsize=9)


# ── chart spec ────────────────────────────────────────────────────────────────
# For each chart: (obs_col, title, exclude_vals, indiv_mode)
# indiv_mode:
#   "per_donor_attr"  — donor has one value for this attribute (Study, Sex, dx, Diagnosis,
#                        Tissue, Tissue_standardized); use mode per donor
#   "donors_with_ct"  — count distinct donors with ≥1 cell of each type (Celltype_markers)

def cell_counts(col, exclude=()):
    s = obs[col].value_counts()
    for ex in exclude:
        s = s.drop(ex, errors="ignore")
    return s

def indiv_counts_attr(col, exclude=()):
    """One value per donor via mode."""
    per = obs.groupby("donor", observed=True)[col].agg(lambda x: x.mode()[0])
    s = per.value_counts()
    for ex in exclude:
        s = s.drop(ex, errors="ignore")
    return s

def indiv_counts_ct(col):
    """Donors with ≥1 cell of each type."""
    s = obs.groupby(col, observed=True)["donor"].nunique().sort_values(ascending=False)
    return s


CHARTS = [
    # (filename_stem, col, title, exclude, indiv_mode, col_key_for_colors)
    ("pie_Celltype_markers", "Celltype_markers", "Cell type",      [],          "donors_with_ct",  "Celltype_markers"),
    ("pie_Study",            "Study",            "Study",          [],          "per_donor_attr",  "Study"),
    ("pie_Sex",              "Sex",              "Sex",            ["unknown"], "per_donor_attr",  "Sex"),
    ("pie_Tissue",           "Tissue",           "Brain Tissue Region", [],    "per_donor_attr",  "Tissue"),
    ("pie_Tissue_standardized", "Tissue_standardized", "Brain region", ["unknown"], "per_donor_attr", "Tissue_standardized"),
    ("pie_Diagnosis",        "Diagnosis",        "Diagnosis",      ["unknown"], "per_donor_attr",  "Diagnosis"),
    ("pie_dx",               "dx",               "Broad diagnosis (dx)", ["unknown"], "per_donor_attr", "dx"),
]

total_cells = len(obs)
total_indiv = obs["donor"].nunique()

for stem, col, title, exclude, indiv_mode, col_key in CHARTS:
    print(f"  {stem} …")

    cc = cell_counts(col, exclude)
    if indiv_mode == "donors_with_ct":
        ic = indiv_counts_ct(col)
        indiv_label = f"n = {ic.sum():,} donor–celltype pairs  |  {total_indiv:,} donors"
        indiv_subtitle = f"n = {total_indiv:,} donors"
    else:
        ic = indiv_counts_attr(col, exclude)
        indiv_subtitle = f"n = {ic.sum():,} individuals"

    cell_subtitle  = f"n = {cc.sum():,} cells"

    colors_c = get_colors(list(cc.index), col_key)
    # keep same colour order for individuals (re-index to cell order)
    ic_aligned = ic.reindex(cc.index).dropna().astype(int)
    ic_aligned = ic_aligned[ic_aligned > 0]
    colors_i   = get_colors(list(ic_aligned.index), col_key)

    fig, axes = plt.subplots(1, 2, figsize=(18, 7))
    fig.patch.set_facecolor("white")

    make_pie(axes[0], cc,          title, cell_subtitle,  colors_c)
    make_pie(axes[1], ic_aligned,  title, indiv_subtitle, colors_i)

    axes[0].set_xlabel("Cells", fontsize=11, labelpad=6)
    axes[1].set_xlabel("Individuals", fontsize=11, labelpad=6)

    fig.tight_layout(rect=[0, 0, 1, 0.97])

    for ext in ("pdf", "png"):
        path = os.path.join(OUT_DIR, f"{stem}.{ext}")
        fig.savefig(path, dpi=180, bbox_inches="tight")
    plt.close(fig)

print(f"\nDone — saved to {OUT_DIR}")
