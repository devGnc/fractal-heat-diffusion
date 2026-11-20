"""Régime établi en fonction de la profondeur du flocon de Koch (même aire).

Usage :
    python3 scripts/profondeur.py                     # CSV de la racine
    python3 scripts/profondeur.py results/v4/regime   # dossier explicite

Lit la ligne de géométrie (« # geometrie=... ») et la dernière ligne de chaque
CSV du dossier (format v3 ou plus récent), affiche un tableau et trace la
conductance de l'enveloppe G = P_ext / (T_air − T_ext), avec T_air la
température moyenne de l'air, en fonction de la profondeur et du périmètre.
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
        return None
    meta = dict(kv.split("=") for kv in first[1:].split())
    if "t_ext" not in meta:
        return None  # ancien format (v2 et avant)
    df = pd.read_csv(path, comment="#")
    last = df.iloc[-1]
    t_ext = float(meta["t_ext"])
    return {
        "geometrie": meta["geometrie"],
        "profondeur": int(meta["profondeur"]),
        "perimetre_m": float(meta["perimetre_m"]),
        "dx_cm": 100 * float(meta["dx"]),
        "T_air": last["T_air"],
        "G": last["P_ext"] / (last["T_air"] - t_ext),
        "bilan_%": 100 * (last["P_ext"] - last["P_int"]) / last["P_int"],
        "jours": last["Time_h"] / 24,
    }


# Les regime_*.csv (make regime, grille fine) sont prioritaires sur les
# séries temporelles (res_*.csv, sessions) du même dossier.
paths = glob.glob(os.path.join(results_dir, "regime_*.csv")) or \
        glob.glob(os.path.join(results_dir, "*.csv"))
rows = [r for p in paths if (r := read_run(p))]
if not rows:
    print(f"Aucun CSV au format v3 ou plus récent dans {results_dir} (lancer make data)")
    sys.exit(1)

table = pd.DataFrame(rows).sort_values(["profondeur"]).set_index("geometrie")
square = table.loc["carre"] if "carre" in table.index else None
if square is not None:
    table["G/G_carre"] = table["G"] / square["G"]
pd.set_option("display.width", 120)
print(table.round(3).to_string())

koch = table[table["profondeur"] >= 0]

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
plt.rcParams["font.size"] = 12

ax1.plot(koch["profondeur"], koch["G"], "o-", color="tab:red", label="Flocon de Koch")
if square is not None:
    ax1.axhline(square["G"], color="black", linestyle="--", label="Carré de même aire")
ax1.set_xlabel("Profondeur de Koch")
ax1.set_ylabel("Conductance de l'enveloppe G (W/(m·K))")
ax1.set_xticks(koch["profondeur"])
ax1.set_title("Chaleur évacuée par degré d'écart")

ax2.plot(koch["perimetre_m"], koch["G"], "o-", color="tab:red", label="Flocon de Koch")
for name, r in koch.iterrows():
    ax2.annotate(f"  n={int(r['profondeur'])}", (r["perimetre_m"], r["G"]))
if square is not None:
    ax2.plot(square["perimetre_m"], square["G"], "ks", label="Carré de même aire")
ax2.set_xlabel("Périmètre intérieur (m)")
ax2.set_ylabel("Conductance de l'enveloppe G (W/(m·K))")
ax2.set_title("G en fonction de la longueur de mur")

for ax in (ax1, ax2):
    ax.grid(True, linestyle="--", alpha=0.7)
    ax.legend()

plt.tight_layout()
plt.show()
