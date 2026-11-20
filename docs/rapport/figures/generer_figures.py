"""Figures du rapport (à lancer depuis la racine du dépôt) :
    python3 docs/rapport/figures/generer_figures.py
Produit koch_profondeurs.png, temoins_formes.png, convergence_test.png, puis les
graphiques des scripts d'analyse (temps.png, profondeur.png, loi_*.png,
temoins.png) à partir de results/v4.
"""
import os
import runpy
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker

OUT = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(OUT))
sys.path.insert(0, os.path.join(REPO, "scripts"))
import formes  # noqa: E402  (mêmes constructions que src/geometry.c)

AREA = 36.0  # m², comme ROOM_AREA


def draw(ax, name, pts, lim=5.6):
    xs = [p[0] for p in pts] + [pts[0][0]]
    ys = [-p[1] for p in pts] + [-pts[0][1]]
    per = formes.perimetre(pts)
    ax.fill(xs, ys, color="#f2c48d", edgecolor="#8a4b12", linewidth=0.6)
    ax.set_title(f"{name}\npérimètre {per:.1f} m".replace(".", ","), fontsize=11)
    ax.set_aspect("equal")
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.axis("off")


# 1. Carré et flocons de profondeur 0 à 4, même aire
fig, axes = plt.subplots(1, 6, figsize=(15, 3.2))
shapes = [("Carré", formes.carre(AREA))] + \
         [(f"Koch {n}", formes.koch(n, AREA)) for n in range(5)]
for ax, (name, pts) in zip(axes, shapes):
    draw(ax, name, pts)
fig.suptitle("Même aire intérieure : 36 m²", fontsize=12)
plt.tight_layout(rect=(0, 0, 1, 0.88))
plt.savefig(os.path.join(OUT, "koch_profondeurs.png"), dpi=170)
plt.close()

# 1 bis. Formes témoins : carrés crénelés (même aire, même périmètre et même
# taille de détail que le flocon de même numéro) et îles de Minkowski
fig, axes = plt.subplots(2, 4, figsize=(12, 6.6))
for n in range(1, 5):
    draw(axes[0][n - 1], f"Carré crénelé {n}", formes.creneau(n, AREA))
for n in range(1, 4):
    draw(axes[1][n - 1], f"Minkowski {n}", formes.minkowski(n, AREA))
axes[1][3].axis("off")
fig.suptitle("Formes témoins, toutes d'aire 36 m²", fontsize=12)
plt.tight_layout(rect=(0, 0, 1, 0.94))
plt.savefig(os.path.join(OUT, "temoins_formes.png"), dpi=170)
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
# 4. Loi de la profondeur optimale et témoins (loi_epaisseur.png,
# loi_seuils.png, temoins.png)
sys.argv = ["scripts/loi.py", OUT]
runpy.run_path("scripts/loi.py", run_name="__main__")
plt.close("all")
print("Figures écrites dans", OUT)
