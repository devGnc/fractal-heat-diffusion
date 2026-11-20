"""Trace la température moyenne en fonction du temps pour le carré et la fractale.

Usage :
    python3 scripts/analyse.py                 # cherche res_*.csv à la racine, sinon results/v0
    python3 scripts/analyse.py results/v0      # dossier explicite
"""
import os
import sys

import matplotlib.pyplot as plt
import pandas as pd

repo_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def find_results_dir():
    if len(sys.argv) > 1:
        return os.path.abspath(sys.argv[1])
    if os.path.exists(os.path.join(repo_dir, "res_carre.csv")):
        return repo_dir
    return os.path.join(repo_dir, "results", "v0")


results_dir = find_results_dir()
path_carre = os.path.join(results_dir, "res_carre.csv")
path_fractal = os.path.join(results_dir, "res_fractal.csv")

try:
    df_carre = pd.read_csv(path_carre)
    df_fractal = pd.read_csv(path_fractal)
except FileNotFoundError as e:
    print(f"Fichier introuvable : {e.filename}")
    sys.exit(1)

print(f"Données lues dans : {results_dir}")

# Les CSV récents ont une colonne P_source (puissance du radiateur) ; pas ceux de v0.
has_power = "P_source" in df_carre and "P_source" in df_fractal
fig, axes = plt.subplots(2 if has_power else 1, 1, figsize=(10, 9 if has_power else 6),
                         sharex=True, squeeze=False)
axes = axes[:, 0]
plt.rcParams["font.size"] = 12

ax = axes[0]
ax.plot(df_carre["Time"], df_carre["Mean"], label="Carré (référence)", color="blue", linewidth=1.5)
ax.plot(df_fractal["Time"], df_fractal["Mean"], label="Fractale (Koch)", color="red", linewidth=1.5)
ax.set_ylabel(r"Température moyenne $T_{moy}$ (°C)")
ax.set_title("Comparaison diffusion thermique : fractale vs carré")

if has_power:
    ax = axes[1]
    ax.plot(df_carre["Time"], df_carre["P_source"], label="Carré (référence)", color="blue", linewidth=1.5)
    ax.plot(df_fractal["Time"], df_fractal["P_source"], label="Fractale (Koch)", color="red", linewidth=1.5)
    ax.set_ylabel("Puissance de la source (u.a.)")
    # Le transitoire initial écrase l'échelle : on cadre sur le régime établi
    tail = pd.concat([df_carre["P_source"].iloc[len(df_carre) // 10:],
                      df_fractal["P_source"].iloc[len(df_fractal) // 10:]])
    ax.set_ylim(0, 1.2 * tail.max())

for ax in axes:
    ax.grid(True, which="both", linestyle="--", alpha=0.7)
    ax.legend()
axes[-1].set_xlabel(r"Temps simulé $t$ (unités de grille)")

plt.tight_layout()
plt.show()
