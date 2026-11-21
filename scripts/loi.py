"""Règle empirique de la profondeur optimale et formes témoins (results/v4/loi/loi.csv).

Usage :
    python3 scripts/loi.py              # tableaux + figures à l'écran
    python3 scripts/loi.py <dossier>    # figures enregistrées dans <dossier>

1. Profondeur optimale en fonction de l'épaisseur e du mur, pour le flocon de
   Koch et l'île de Minkowski.
2. Seuils : épaisseur e_c(n) à laquelle la profondeur n cesse de battre la
   profondeur n − 1 (G_n = G_{n−1}), comparée à la taille ℓ_n des détails
   ajoutés à la profondeur n, pour trois aires.
3. Témoins : à e = 20 cm, flocon, carré crénelé de même périmètre, Minkowski.
4. Ce qui compte : G en fonction de la surface extérieure du mur.
"""
import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import formes  # noqa: E402

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
out_dir = os.path.abspath(sys.argv[1]) if len(sys.argv) > 1 else None

df = pd.read_csv(os.path.join(REPO, "results", "v4", "loi", "loi.csv"))
df["detail_m"] = [formes.detail(f, n, a) if f != "carre" else np.nan
                  for f, n, a in zip(df["famille"], df["profondeur"], df["aire_m2"])]
sq = df[df["famille"] == "carre"].set_index(["aire_m2", "epaisseur_m"])["G"]
df["G_rel"] = df["G"] / sq.reindex(list(zip(df["aire_m2"], df["epaisseur_m"]))).values
pd.set_option("display.width", 160)


def save(name):
    if out_dir:
        plt.savefig(os.path.join(out_dir, name), dpi=160)
        plt.close()


def serie(famille, aire):
    """Tableau G[e, profondeur] ; le carré sert de profondeur 0 à Minkowski."""
    d = df[(df["aire_m2"] == aire) & (df["famille"] == famille)]
    t = d.pivot_table(index="epaisseur_m", columns="profondeur", values="G")
    if famille == "minkowski":
        t[0] = sq.xs(aire, level="aire_m2")
        t = t[sorted(t.columns)]
    return t


def seuil(t, n):
    """Première épaisseur où G_n − G_{n−1} change de signe (interpolation
    linéaire), ou None."""
    if n not in t or n - 1 not in t:
        return None
    diff = (t[n] - t[n - 1]).dropna()
    e, v = diff.index.values, diff.values
    for i in range(1, len(v)):
        if v[i - 1] > 0 >= v[i]:
            return e[i - 1] + (e[i] - e[i - 1]) * v[i - 1] / (v[i - 1] - v[i])
    return None


# 1. Profondeur optimale --------------------------------------------------
print("Profondeur optimale (36 m²) en fonction de l'épaisseur du mur :")
opt = pd.DataFrame({f: serie(f, 36.0).idxmax(axis=1) for f in ("koch", "minkowski")})
print(opt.T.to_string(), "\n")

fig, axes = plt.subplots(1, 2, figsize=(13, 4.8), sharey=True)
colors = plt.cm.viridis(np.linspace(0, 0.9, 5))
for ax, famille, titre in ((axes[0], "koch", "Flocon de Koch"),
                           (axes[1], "minkowski", "Île de Minkowski")):
    t = serie(famille, 36.0)
    for n in t.columns:
        if n == 0 and famille == "minkowski":
            continue
        c = colors[n]
        ax.plot(t.index * 100, t[n] / sq.xs(36.0, level="aire_m2").reindex(t.index),
                "o-", ms=3, color=c, label=f"profondeur {n}")
        ell = formes.detail(famille, n, 36.0)
        if n >= 1 and ell <= 0.5:
            ax.axvline(ell * 100, color=c, ls=":", lw=1.2)
            ax.annotate(f"ℓ{n} = {ell * 100:.0f} cm".replace(".", ","),
                        (ell * 100, 0.98), xycoords=("data", "axes fraction"),
                        rotation=90, va="top", ha="right", fontsize=8, color=c)
    ax.axhline(1, color="k", ls="--", lw=1)
    ax.set_xscale("log")
    ax.set_xlim(4.5, 50)
    ax.set_xticks([5, 10, 20, 30, 45])
    ax.set_xticklabels(["5", "10", "20", "30", "45"])
    ax.set_xlabel("Épaisseur du mur e (cm)")
    ax.set_title(f"{titre}, 36 m²")
    ax.grid(True, which="both", ls="--", alpha=0.4)
    ax.legend(fontsize=8)
axes[0].set_ylabel("G / G du carré")
plt.tight_layout()
save("loi_epaisseur.png")

# 2. Seuils ----------------------------------------------------------------
rows = []
for famille in ("koch", "minkowski"):
    for aire in sorted(df["aire_m2"].unique()):
        t = serie(famille, aire)
        for n in range(1, int(t.columns.max()) + 1):
            ec = seuil(t, n)
            if ec is not None:
                ell = formes.detail(famille, n, aire)
                rows.append({"famille": famille, "aire_m2": aire, "n": n,
                             "detail_cm": 100 * ell, "seuil_cm": 100 * ec,
                             "seuil/detail": ec / ell})
seuils = pd.DataFrame(rows)
print("Seuils : la profondeur n cesse de battre n − 1 quand e > e_c(n)")
print(seuils.round(3).to_string(index=False))
for famille, g in seuils.groupby("famille"):
    print(f"  {famille} : e_c/ℓ = {g['seuil/detail'].mean():.2f} "
          f"(de {g['seuil/detail'].min():.2f} à {g['seuil/detail'].max():.2f})")
print()

fig, ax = plt.subplots(figsize=(6.2, 5.2))
marks = {16.0: "v", 36.0: "o", 64.0: "^"}
for (famille, aire), g in seuils.groupby(["famille", "aire_m2"]):
    ax.loglog(g["detail_cm"], g["seuil_cm"], marks[aire], ms=7,
              color="tab:red" if famille == "koch" else "tab:blue",
              label=f"{'Koch' if famille == 'koch' else 'Minkowski'}, {aire:g} m²")
x = np.array([4, 150])
for famille, g in seuils.groupby("famille"):
    c = g["seuil/detail"].mean()
    ax.loglog(x, c * x, "--", lw=1, color="tab:red" if famille == "koch" else "tab:blue",
              label=f"e_c = {c:.2f} ℓ".replace(".", ","))
ax.loglog(x, x, "k:", lw=1, label="e_c = ℓ")
ax.set_xlabel("Taille ℓ des détails ajoutés (cm)")
ax.set_ylabel("Épaisseur seuil e_c (cm)")
ax.set_title("Un niveau de détail sert tant que le mur\nest plus fin que lui")
ax.grid(True, which="both", ls="--", alpha=0.4)
ax.legend(fontsize=8)
plt.tight_layout()
save("loi_seuils.png")

# 3. Témoins à e = 20 cm ---------------------------------------------------
ref = df[(df["aire_m2"] == 36.0) & (df["epaisseur_m"] == 0.20)].copy()
print("Témoins, 36 m², mur de 20 cm :")
print(ref.sort_values(["famille", "profondeur"])[
    ["geometrie", "perimetre_m", "detail_m", "surface_ext_m", "T_air", "G", "G_rel"]
].round(3).to_string(index=False), "\n")

styles = {"carre": ("s", "Carré"), "koch": ("o-", "Flocon de Koch"),
          "creneau": ("D--", "Carré crénelé (même périmètre que Koch n)"),
          "minkowski": ("^-", "Île de Minkowski")}
pal = {"carre": "k", "koch": "tab:red", "creneau": "tab:gray", "minkowski": "tab:blue"}
fig, axes = plt.subplots(1, 2, figsize=(13, 4.8))
for famille, (m, label) in styles.items():
    g = ref[ref["famille"] == famille].sort_values("profondeur")
    for ax, xcol in zip(axes, ("perimetre_m", "surface_ext_m")):
        ax.plot(g[xcol], g["G"], m, color=pal[famille], label=label)
        if famille != "carre":
            for _, r in g.iterrows():
                ax.annotate(f" {int(r['profondeur'])}", (r[xcol], r["G"]),
                            fontsize=8, color=pal[famille])
axes[0].set_xlabel("Périmètre intérieur (m)")
axes[1].set_xlabel("Surface extérieure du mur (m, par mètre de hauteur)")
for ax in axes:
    ax.set_ylabel("Conductance G (W/(m·K))")
    ax.grid(True, ls="--", alpha=0.4)
axes[0].set_title("Même périmètre, résultats très différents")
axes[1].set_title("La surface extérieure du mur ordonne tout")
axes[0].legend(fontsize=8)
plt.tight_layout()
save("temoins.png")

# 4. Toutes les formes, toutes les épaisseurs : G et surface extérieure ----
d36 = df[df["aire_m2"] == 36.0]
corr = np.corrcoef(np.log(d36["surface_ext_m"]), np.log(d36["G"]))[0, 1]
corr_p = np.corrcoef(np.log(d36["perimetre_m"]), np.log(d36["G"]))[0, 1]
print(f"Corrélation (log) de G, à 36 m² et toutes épaisseurs confondues, avec :"
      f" la surface extérieure {corr:.3f}, le périmètre intérieur {corr_p:.3f}")

if not out_dir:
    plt.show()
