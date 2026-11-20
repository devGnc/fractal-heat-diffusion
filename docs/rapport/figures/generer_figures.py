"""Figures du rapport (à lancer depuis la racine du dépôt) :
    python3 docs/rapport/figures/generer_figures.py
Produit koch_profondeurs.png, convergence_test.png, puis les graphiques des
scripts d'analyse (temps.png, profondeur.png) à partir de results/v4.
"""
import math
import os
import runpy
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker

OUT = os.path.dirname(os.path.abspath(__file__))
AREA = 36.0  # m², comme ROOM_AREA


def koch_side(a, b, depth, out):
    """Même construction que koch_side() dans src/geometry.c."""
    if depth == 0:
        out.append(a)
        return
    dx, dy = (b[0] - a[0]) / 3, (b[1] - a[1]) / 3
    p3 = (a[0] + dx, a[1] + dy)
    p4 = (a[0] + 2 * dx, a[1] + 2 * dy)
    ang = -math.pi / 3
    tip = (p3[0] + dx * math.cos(ang) - dy * math.sin(ang),
           p3[1] + dx * math.sin(ang) + dy * math.cos(ang))
    for u, v in ((a, p3), (p3, tip), (tip, p4), (p4, b)):
        koch_side(u, v, depth - 1, out)


def snowflake(depth):
    f = 8 / 5 - 3 / 5 * (4 / 9) ** depth
    r = math.sqrt(AREA / (f * 3 * math.sqrt(3) / 4))
    s = r * math.sqrt(3)
    p1, p2, p3 = (-s / 2, r / 2), (s / 2, r / 2), (0, -r)
    pts = []
    for u, v in ((p1, p3), (p3, p2), (p2, p1)):
        koch_side(u, v, depth, pts)
    return pts


# 1. Carré et flocons de profondeur 0 à 4, même aire
fig, axes = plt.subplots(1, 6, figsize=(15, 3.2))
h = math.sqrt(AREA) / 2
shapes = [("Carré", [(-h, -h), (h, -h), (h, h), (-h, h)])] + \
         [(f"Koch {n}", snowflake(n)) for n in range(5)]
for ax, (name, pts) in zip(axes, shapes):
    xs = [p[0] for p in pts] + [pts[0][0]]
    ys = [-p[1] for p in pts] + [-pts[0][1]]
    per = sum(math.dist(pts[i], pts[(i + 1) % len(pts)]) for i in range(len(pts)))
    ax.fill(xs, ys, color="#f2c48d", edgecolor="#8a4b12", linewidth=0.8)
    ax.set_title(f"{name}\npérimètre {per:.1f} m".replace(".", ","), fontsize=11)
    ax.set_aspect("equal")
    ax.set_xlim(-5.6, 5.6)
    ax.set_ylim(-5.6, 5.6)
    ax.axis("off")
fig.suptitle("Même aire intérieure : 36 m²", fontsize=12)
plt.tight_layout(rect=(0, 0, 1, 0.88))
plt.savefig(os.path.join(OUT, "koch_profondeurs.png"), dpi=170)
plt.close()

# 2. Convergence du test de la pièce circulaire (valeurs de make test)
dx = [5.0, 2.5, 1.25]
cases = {
    "Sans résistances de surface (corrigé)": [0.548, 0.276, 0.159],
    "Surfaces de la norme, corrigées": [1.512, 0.847, 0.453],
    "Surfaces de la norme, sans correction": [13.505, 13.778, 13.731],
}
fig, ax = plt.subplots(figsize=(6.4, 4.4))
for (label, err), m in zip(cases.items(), ["o-", "s-", "^--"]):
    ax.loglog(dx, err, m, label=label)
ax.loglog(dx, [0.3 * d for d in dx], "k:", label="pente 1 (ordre 1)")
ax.set_xlabel("Pas de la grille Δx (cm)")
ax.set_ylabel("|erreur sur G| (%)")
ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
ax.set_xticks(dx)
ax.set_xticklabels(["5", "2,5", "1,25"])
ax.grid(True, which="both", linestyle="--", alpha=0.5)
ax.legend(fontsize=9)
ax.set_title("Test de la pièce circulaire")
plt.tight_layout()
plt.savefig(os.path.join(OUT, "convergence_test.png"), dpi=170)
plt.close()

# 3. Graphiques des scripts d'analyse
for script, arg, name in (("analyse.py", "results/v4/temps", "temps.png"),
                          ("profondeur.py", "results/v4/regime", "profondeur.png")):
    plt.show = lambda n=name: plt.savefig(os.path.join(OUT, n), dpi=150)
    sys.argv = ["scripts/" + script, arg]
    runpy.run_path("scripts/" + script, run_name="__main__")
    plt.close("all")
print("Figures écrites dans", OUT)
