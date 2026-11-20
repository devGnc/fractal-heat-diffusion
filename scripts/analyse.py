"""Trace l'évolution temporelle de chaque simulation.

Usage :
    python3 scripts/analyse.py                 # res_*.csv à la racine, sinon results/v3
    python3 scripts/analyse.py results/v2      # dossier explicite

Chaque fichier res_<geometrie>.csv du dossier donne une courbe. Format v3 :
température de l'air et puissance évacuée par les murs en fonction du temps
(heures). Les anciens formats (v0 à v2, unités de grille) restent lisibles.
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
    return os.path.join(repo_dir, "results", "v3")


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

v3 = all("T_air" in df for df in runs.values())
if v3:
    time_col, temp_col, power_col = "Time_h", "T_air", "P_ext"
    temp_label = r"Température de l'air $T_{air}$ (°C)"
    power_label = "Puissance évacuée par les murs (W/m)"
    time_label = "Temps simulé (h)"
else:
    time_col, temp_col, power_col = "Time", "Mean", "P_source"
    temp_label = r"Température moyenne de l'air $T_{moy}$ (°C)"
    power_label = "Puissance de la source (u.a.)"
    time_label = r"Temps simulé $t$ (unités de grille)"
has_power = all(power_col in df for df in runs.values())  # absent en v0

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
    axes[0].plot(df[time_col], df[temp_col], label=label(name), **style[name])
axes[0].set_ylabel(temp_label)
axes[0].set_title("Évolution vers le régime établi")

if has_power:
    for name in names:
        df = runs[name]
        axes[1].plot(df[time_col], df[power_col], label=label(name), **style[name])
    axes[1].set_ylabel(power_label)
    if not v3:
        # Le transitoire initial écrase l'échelle : on cadre sur le régime établi
        tail = pd.concat([df[power_col].iloc[len(df) // 10:] for df in runs.values()])
        axes[1].set_ylim(0, 1.2 * tail.max())

for ax in axes:
    ax.grid(True, which="both", linestyle="--", alpha=0.7)
    ax.legend()
axes[-1].set_xlabel(time_label)

plt.tight_layout()
plt.show()
