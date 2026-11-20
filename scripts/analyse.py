"""Trace l'évolution temporelle de chaque simulation (T moyenne, puissance source).

Usage :
    python3 scripts/analyse.py                 # res_*.csv à la racine, sinon results/v0
    python3 scripts/analyse.py results/v2      # dossier explicite

Chaque fichier res_<geometrie>.csv du dossier donne une courbe.
"""
import glob
import os
import sys

import matplotlib.pyplot as plt
import pandas as pd

repo_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

LABELS = {"carre": "Carré iso-aire", "fractal": "Fractale (Koch)"}


def find_results_dir():
    if len(sys.argv) > 1:
        return os.path.abspath(sys.argv[1])
    if glob.glob(os.path.join(repo_dir, "res_*.csv")):
        return repo_dir
    return os.path.join(repo_dir, "results", "v0")


def sort_key(name):
    # Carré d'abord, puis les profondeurs de Koch dans l'ordre
    return (name != "carre", name)


def label(name):
    if name.startswith("koch"):
        return f"Koch profondeur {name[4:]}"
    return LABELS.get(name, name)


results_dir = find_results_dir()
paths = glob.glob(os.path.join(results_dir, "res_*.csv"))
if not paths:
    print(f"Aucun fichier res_*.csv dans {results_dir}")
    sys.exit(1)

runs = {}
for path in paths:
    name = os.path.basename(path)[len("res_"):-len(".csv")]
    runs[name] = pd.read_csv(path, comment="#")
names = sorted(runs, key=sort_key)
print(f"Données lues dans : {results_dir} ({', '.join(names)})")

# Les anciens CSV (v0) n'ont pas de colonne P_source
has_power = all("P_source" in df for df in runs.values())
fig, axes = plt.subplots(2 if has_power else 1, 1, figsize=(10, 9 if has_power else 6),
                         sharex=True, squeeze=False)
axes = axes[:, 0]
plt.rcParams["font.size"] = 12

colors = plt.cm.viridis([i / max(1, len(names) - 1) for i in range(len(names))])
style = {name: dict(color="black" if name == "carre" else colors[i],
                    linestyle="--" if name == "carre" else "-", linewidth=1.5)
         for i, name in enumerate(names)}

for name in names:
    df = runs[name]
    axes[0].plot(df["Time"], df["Mean"], label=label(name), **style[name])
axes[0].set_ylabel(r"Température moyenne de l'air $T_{moy}$ (°C)")
axes[0].set_title("Diffusion thermique : évolution vers le régime établi")

if has_power:
    for name in names:
        df = runs[name]
        axes[1].plot(df["Time"], df["P_source"], label=label(name), **style[name])
    axes[1].set_ylabel("Puissance de la source (u.a.)")
    # Le transitoire initial écrase l'échelle : on cadre sur le régime établi
    tail = pd.concat([df["P_source"].iloc[len(df) // 10:] for df in runs.values()])
    axes[1].set_ylim(0, 1.2 * tail.max())

for ax in axes:
    ax.grid(True, which="both", linestyle="--", alpha=0.7)
    ax.legend()
axes[-1].set_xlabel(r"Temps simulé $t$ (unités de grille)")

plt.tight_layout()
plt.show()
