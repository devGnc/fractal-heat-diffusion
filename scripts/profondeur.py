"""Régime établi en fonction de la profondeur du flocon de Koch (même aire).

Usage :
    python3 scripts/profondeur.py              # res_*.csv à la racine
    python3 scripts/profondeur.py results/v2   # dossier explicite

Lit la ligne de géométrie (« # geometrie=... ») et la dernière ligne de chaque
res_koch<n>.csv et res_carre.csv, affiche un tableau et trace la puissance de
la source en fonction de la profondeur. Comparer les profondeurs de Koch entre
elles neutralise le biais de grille des murs obliques (mêmes orientations
0°, 60°, 120° à toutes les profondeurs).
"""
import glob
import os
import sys

import matplotlib.pyplot as plt
import pandas as pd

repo_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
results_dir = os.path.abspath(sys.argv[1]) if len(sys.argv) > 1 else repo_dir


def read_run(path):
    with open(path) as f:
        first = f.readline()
    if not first.startswith("#"):
        return None  # ancien format, sans géométrie
    meta = dict(kv.split("=") for kv in first[1:].split())
    df = pd.read_csv(path, comment="#")
    last = df.iloc[-1]
    return {
        "geometrie": meta["geometrie"],
        "profondeur": int(meta["profondeur"]),
        "aire": int(meta["aire"]),
        "faces_air_mur": int(meta["faces_air_mur"]),
        "T_moy": last["Mean"],
        "T_min": last["Min"],
        "P_source": last["P_source"],
        "bilan_%": 100 * (last["P_ext"] + last["P_cool"] - last["P_source"]) / last["P_source"],
        "t_fin": last["Time"],
    }


rows = [r for p in glob.glob(os.path.join(results_dir, "res_*.csv")) if (r := read_run(p))]
if not rows:
    print(f"Aucun CSV au nouveau format dans {results_dir} (lancer make data)")
    sys.exit(1)

table = pd.DataFrame(rows).sort_values(["profondeur"]).set_index("geometrie")
pd.set_option("display.width", 120)
print(table.round(3).to_string())

koch = table[table["profondeur"] >= 0]
square = table.loc["carre"] if "carre" in table.index else None

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
plt.rcParams["font.size"] = 12

ax1.plot(koch["profondeur"], koch["P_source"], "o-", color="tab:red", label="Flocon de Koch")
if square is not None:
    ax1.axhline(square["P_source"], color="black", linestyle="--", label="Carré iso-aire")
ax1.set_xlabel("Profondeur de Koch")
ax1.set_ylabel("Puissance de la source en régime établi (u.a.)")
ax1.set_xticks(koch["profondeur"])
ax1.set_title("Puissance nécessaire pour maintenir la source")

ax2.plot(koch["faces_air_mur"], koch["P_source"], "o-", color="tab:red", label="Flocon de Koch")
for name, r in koch.iterrows():
    ax2.annotate(f"  n={int(r['profondeur'])}", (r["faces_air_mur"], r["P_source"]))
if square is not None:
    ax2.plot(square["faces_air_mur"], square["P_source"], "ks", label="Carré iso-aire")
ax2.set_xlabel("Surface d'échange : nombre de faces air/mur")
ax2.set_ylabel("Puissance de la source (u.a.)")
ax2.set_title("Puissance en fonction de la surface d'échange")

for ax in (ax1, ax2):
    ax.grid(True, linestyle="--", alpha=0.7)
    ax.legend()

plt.tight_layout()
plt.show()
