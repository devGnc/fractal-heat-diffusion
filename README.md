# TIPE : diffusion thermique et géométries fractales

TIPE MPI de Devrim Genc et Mathis Pagnier. Ce dépôt reprend la simulation de
l'ancien dépôt GitLab du lycée (`tests/diffusion`).

## Problématique

> Comment utiliser les fractales pour optimiser la diffusion thermique en architecture ?

On compare une pièce dont le contour est un **flocon de Koch** à une **pièce carrée
de même aire**, chauffées par une source centrale à 100 °C, avec l'extérieur maintenu
à 20 °C. À aire égale, la fractale a un périmètre bien plus long : elle échange donc
davantage de chaleur avec l'extérieur.

Mots-clés : fractales, suites, récursivité, programmation en C, équation de la
chaleur, différences finies, complexité.

## Compilation et exécution

```bash
brew install sdl2 pkg-config      # macOS
make                              # compile ./thermalfractalsim
make run                          # fenêtre interactive
make data                         # simulations complètes sans fenêtre (≈ 25 s)
make plot                         # trace les courbes (pandas + matplotlib)
```

Mode sans fenêtre, pour obtenir des données exploitables (enregistrement dès
l'itération 0, une ligne toutes les 400 itérations) :

```bash
./thermalfractalsim --headless carre                 # 300 000 itérations → res_carre.csv
./thermalfractalsim --headless fractal 100000 f.csv  # nombre d'itérations et fichier choisis
make data ITER=200000
```

| Touche | Action |
|---|---|
| Espace | Basculer carré ↔ fractale |
| S | Démarrer / arrêter l'enregistrement CSV (`res_carre.csv` ou `res_fractal.csv`) ; Espace et R l'arrêtent aussi |
| R | Réinitialiser |
| Q | Quitter |

Colonnes du CSV : `Iter,Time,Mean,Std,Min,P_source`.

- `Mean`, `Std`, `Min` : température de l'**air intérieur**, murs et source exclus.
- `P_source` : puissance fournie par le radiateur pour rester à 100 °C, c'est-à-dire
  l'énergie réinjectée dans la source à chaque pas divisée par Δt (°C·px par unité
  de temps, capacité thermique 1 par pixel). En régime établi, c'est la puissance
  perdue par la pièce.

Les CSV de `results/v0` ont l'ancien format (`Iter,Time,Mean,Std,Max`, murs et
source inclus dans la moyenne). Les CSV produits à la racine sont ignorés par git :
il faut les déplacer dans `results/<version>/` pour les garder.

## Organisation

```
src/
  thermal_diffusion.h   paramètres physiques et structure Sim
  geometry.c            Koch récursif, Bresenham, flood fill, carré iso-aire
  physics.c             schéma explicite (FTCS), CFL, statistiques
  render.c              carte de couleurs SDL
  main.c                boucle SDL, clavier, enregistrement CSV
scripts/analyse.py      tracé T_moy(t)
results/v0/             premiers résultats (ancien dépôt, régime non établi)
results/v1/             mode headless, stats sur l'air, schéma encore non conservatif
docs/                   export GeoGebra
```

## Modèle actuel (v0)

- Grille 300×300, Δx = 1, Δt = 0,05, α_air = 0,20, α_mur = 0,02
- Extérieur à 20 °C (Dirichlet), source disque de rayon 15 à 100 °C
- Pertes de Newton faibles (`COOLING_FACTOR`)
- Koch de profondeur 2 ; carré de même aire (approximativement)

`make data` (300 000 itérations ; le régime est établi dès 200 000 environ) :

| | Aire intérieure | Pixels de mur | T moyenne (air, hors source) | Puissance source |
|---|---|---|---|---|
| Fractale | 13 389 px | 698 | 37,99 °C | 77,75 |
| Carré | 12 769 px | 456 | 41,68 °C | 73,12 |

## Limites connues et pistes

1. **Objectif à préciser** : la fractale isole moins bien mais dissipe mieux
   (radiateurs, dissipateurs, façades qui évacuent la chaleur).
2. **Iso-aire approximative** : le carré est 4,6 % plus petit que la fractale (`sqrtf` tronqué).
3. ~~Murs inclus dans la moyenne~~ : corrigé, la moyenne porte sur l'air hors source.
4. **Schéma non conservatif** quand α varie : utiliser ∇·(α∇T) avec la moyenne
   harmonique de α aux interfaces. Conséquence : en régime établi, T vérifie ΔT ≈ 0
   quel que soit α, donc **les murs n'isolent pas** (T_min ≈ 20,0 °C contre le mur).
5. ~~Régime non établi~~ dans `results/v0` : corrigé par le mode `--headless`.
6. **CFL très large** (0,01 pour un maximum de 0,25) : Δt peut être environ 20 fois plus grand.
7. **Unités de grille** : adimensionner ou fixer des Δx et Δt physiques.
8. **Mesure plus parlante** : `P_source` en régime stationnaire en fonction de la
   profondeur fractale (0, 1, 2, 3).
9. **Biais de grille** : avec le stencil à 5 points, un mur oblique en escalier
   échange sur sa longueur |Δx|+|Δy|. La fractale compte 934 faces air/mur pour un
   périmètre géométrique d'environ 792 (+18 %), le carré 452 (exact). Ce biais ne
   diminue pas en raffinant la grille ; comparer les profondeurs de Koch entre elles
   l'annule (mêmes orientations 0°, 60°, 120°).
