"""Balayage en épaisseur de mur et en aire : régime établi sur la grille fine.

Usage (depuis la racine, après make build/thermalfractalsim_fin) :
    python3 scripts/balayage.py            # complète results/v4/loi/loi.csv
    python3 scripts/balayage.py -j 4       # 4 calculs en parallèle (défaut : cœurs − 1)

Chaque calcul est un régime établi (--regime) d'une géométrie, pour une épaisseur
de mur e et une aire A. Le résultat est ajouté à results/v4/loi/loi.csv (une
ligne par calcul) ; les calculs déjà présents sont sautés, si bien qu'un
balayage interrompu reprend là où il s'était arrêté. Les géométries qui ne
tiennent pas dans la grille sont signalées et ignorées.
"""
import argparse
import csv
import os
import subprocess
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor, as_completed

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BIN = os.path.join(REPO, "build", "thermalfractalsim_fin")
OUT = os.path.join(REPO, "results", "v4", "loi", "loi.csv")

# Épaisseurs de mur réalistes, de 5 à 45 cm, resserrées là où se trouvent les
# seuils (de quelques cm à une vingtaine de cm à 36 m²).
EPAISSEURS = [0.05, 0.06, 0.07, 0.08, 0.09, 0.10, 0.11, 0.12, 0.14, 0.16, 0.18,
              0.20, 0.22, 0.24, 0.26, 0.28, 0.30, 0.32, 0.35, 0.40, 0.45]
KOCH = ["carre"] + [f"koch{n}" for n in range(5)]
MINKOWSKI = [f"minkowski{n}" for n in range(1, 4)]
CRENEAU = [f"creneau{n}" for n in range(1, 5)]

# (aire, géométries) : toutes les formes à 36 m². D'autres aires peuvent être
# ajoutées ici (par exemple (16.0, KOCH + MINKOWSKI)) pour vérifier que le
# seuil suit la taille des détails ; loi.py les prend en compte.
PLAN = [(36.0, KOCH + CRENEAU + MINKOWSKI)]

COLONNES = ["aire_m2", "epaisseur_m", "geometrie", "famille", "profondeur",
            "perimetre_m", "surface_int_m", "surface_ext_m", "dx_m",
            "T_air", "T_air_max", "T_mur", "P_ext", "G", "bilan_pct"]


def famille(g):
    for f in ("koch", "creneau", "minkowski"):
        if g.startswith(f):
            return f
    return "carre"


def calcul(aire, e, g):
    """Lance un régime établi ; renvoie la ligne de résultats, ou None."""
    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, "r.csv")
        p = subprocess.run([BIN, "--regime", g, path, "--epaisseur", str(e),
                            "--aire", str(aire)], capture_output=True, text=True)
        if p.returncode != 0:
            return None, p.stderr.strip().splitlines()[-1] if p.stderr else "échec"
        with open(path) as f:
            meta = dict(kv.split("=") for kv in f.readline()[1:].split())
            next(f)
            last = list(csv.reader(f))[-1]
    T_air, T_max, T_mur, P_int, P_mur, P_ext = map(float, last[2:8])
    t_ext = float(meta["t_ext"])
    return {
        "aire_m2": aire, "epaisseur_m": e, "geometrie": g, "famille": famille(g),
        "profondeur": int(meta["profondeur"]),
        "perimetre_m": float(meta["perimetre_m"]),
        "surface_int_m": float(meta["surface_int_m"]),
        "surface_ext_m": float(meta["surface_ext_m"]),
        "dx_m": float(meta["dx"]),
        "T_air": T_air, "T_air_max": T_max, "T_mur": T_mur, "P_ext": P_ext,
        "G": P_ext / (T_air - t_ext),
        "bilan_pct": 100 * (P_ext - P_int) / P_int,
    }, None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-j", type=int, default=max(1, (os.cpu_count() or 2) - 1))
    args = ap.parse_args()
    if not os.path.exists(BIN):
        sys.exit(f"{BIN} absent : lancer make build/thermalfractalsim_fin")

    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    done = set()
    if os.path.exists(OUT):
        with open(OUT) as f:
            for r in csv.DictReader(f):
                try:  # une ligne tronquée (balayage tué en pleine écriture) est refaite
                    float(r["bilan_pct"])
                    done.add((float(r["aire_m2"]), float(r["epaisseur_m"]), r["geometrie"]))
                except (TypeError, ValueError):
                    pass
    todo = [(a, e, g) for a, geoms in PLAN for e in EPAISSEURS for g in geoms
            if (a, e, g) not in done]
    print(f"{len(done)} calculs déjà faits, {len(todo)} à faire ({args.j} en parallèle)")

    new = not os.path.exists(OUT)
    with open(OUT, "a", newline="") as f, ThreadPoolExecutor(args.j) as pool:
        w = csv.DictWriter(f, COLONNES)
        if new:
            w.writeheader()
        futures = {pool.submit(calcul, *t): t for t in todo}
        for k, fut in enumerate(as_completed(futures), 1):
            a, e, g = futures[fut]
            row, err = fut.result()
            if row is None:
                print(f"[{k}/{len(todo)}] A={a:g} e={e:g} {g} : ignoré ({err})")
                continue
            w.writerow({c: (f"{row[c]:.6g}" if isinstance(row[c], float) else row[c])
                        for c in COLONNES})
            f.flush()
            print(f"[{k}/{len(todo)}] A={a:g} e={e:g} {g} : G = {row['G']:.2f}")


if __name__ == "__main__":
    main()
