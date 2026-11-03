#!/usr/bin/env python3
"""Swatch sheet of the paper-wide colour scheme in src/viz/style.py, for BioRender (Fig 1)
and for review. Reads every colour from style.py, so the sheet cannot drift from the scheme.
Output: figures/palette_reference.png / .pdf (7.2 x 5 in, 600 dpi PNG + vector PDF).
"""
import sys
sys.path.insert(0, "src")
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from viz.style import (apply_style, save_fig, METHOD_COLORS, MOLECULE_COLORS, REGION_COLORS,
                       PROTOCOL_COLORS, CMAP_SEQ)

apply_style()

SECTIONS = [("Methods", METHOD_COLORS), ("Molecules", MOLECULE_COLORS),
            ("MHC regions", REGION_COLORS), ("Protocols", PROTOCOL_COLORS)]
W, H, GAP = 0.121, 0.071, 0.0135   # swatch width/height and gap, figure fractions

fig = plt.figure(figsize=(7.2, 5.0))
y = 0.985
for title, colors in SECTIONS:
    fig.text(0.017, y, title, fontsize=8, fontweight="bold", va="top")
    top = y - 0.042
    for i, (name, hexc) in enumerate(colors.items()):
        x0 = 0.03 + i * (W + GAP)
        fig.patches.append(Rectangle((x0, top - H), W, H, transform=fig.transFigure,
                                     facecolor=hexc, edgecolor="0.6", lw=0.5))
        fig.text(x0 + W / 2, top - H - 0.012, name, fontsize=7, ha="center", va="top")
        fig.text(x0 + W / 2, top - H - 0.035, hexc.upper(), fontsize=6.5, ha="center",
                 va="top", color="0.4")
    y -= 0.216
fig.text(0.017, y, "Sequential heatmap (all heatmaps + PyMOL footprint)", fontsize=8,
         fontweight="bold", va="top")
ax = fig.add_axes([0.017, y - 0.042 - H * 0.8, 0.966, H * 0.8])
ax.imshow(np.linspace(0, 1, 512)[None, :], aspect="auto", cmap=CMAP_SEQ)
ax.set_xticks([]); ax.set_yticks([])

save_fig(fig, "figures/palette_reference")
plt.close(fig)
for title, colors in SECTIONS:
    print(f"{title}: " + ", ".join(f"{k} {v}" for k, v in colors.items()))
