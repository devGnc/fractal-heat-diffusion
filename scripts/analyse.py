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

plt.figure(figsize=(10, 6))
plt.rcParams["font.size"] = 12

plt.plot(df_carre["Time"], df_carre["Mean"], label="Carré (référence)", color="blue", linewidth=1.5)
plt.plot(df_fractal["Time"], df_fractal["Mean"], label="Fractale (Koch)", color="red", linewidth=1.5)

plt.xlabel(r"Temps simulé $t$ (unités de grille)")
plt.ylabel(r"Température moyenne $T_{moy}$ (°C)")
plt.title("Comparaison diffusion thermique : fractale vs carré")
plt.grid(True, which="both", linestyle="--", alpha=0.7)
plt.legend()

plt.tight_layout()
plt.show()
