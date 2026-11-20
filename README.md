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
make run                          # lance la simulation
make plot                         # trace les courbes (pandas + matplotlib)
```

| Touche | Action |
|---|---|
| Espace | Basculer carré ↔ fractale |
| S | Démarrer / arrêter l'enregistrement CSV (`res_carre.csv` ou `res_fractal.csv`) |
| R | Réinitialiser |
| Q | Quitter |

Le CSV contient les colonnes `Iter,Time,Mean,Std,Max`. Les CSV produits à la racine
sont ignorés par git : il faut les déplacer dans `results/<version>/` pour les garder.

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
docs/                   export GeoGebra
```

## Modèle actuel (v0)

- Grille 300×300, Δx = 1, Δt = 0,05, α_air = 0,20, α_mur = 0,02
- Extérieur à 20 °C (Dirichlet), source disque de rayon 15 à 100 °C
- Pertes de Newton faibles (`COOLING_FACTOR`)
- Koch de profondeur 2 ; carré de même aire (approximativement)

Vérification sans fenêtre graphique, après 300 000 itérations (régime établi) :

| | Aire intérieure | Pixels de mur | T moyenne intérieure |
|---|---|---|---|
| Fractale | 13 389 px | 698 | 41,3 °C |
| Carré | 12 769 px | 456 | 44,9 °C |

## Limites connues et pistes

1. **Objectif à préciser** : la fractale isole moins bien mais dissipe mieux
   (radiateurs, dissipateurs, façades qui évacuent la chaleur).
2. **Iso-aire approximative** : le carré est 4,6 % plus petit que la fractale (`sqrtf` tronqué).
3. **Murs inclus dans la moyenne** (`compute_stats`), ce qui pénalise la fractale.
4. **Schéma non conservatif** quand α varie : utiliser ∇·(α∇T) avec la moyenne
   harmonique de α aux interfaces.
5. **Régime non établi** dans `results/v0` ; les deux enregistrements démarrent à
   des instants différents.
6. **CFL très large** (0,01 pour un maximum de 0,25) : Δt peut être environ 20 fois plus grand.
7. **Unités de grille** : adimensionner ou fixer des Δx et Δt physiques.
8. **Mesure plus parlante** : flux de chaleur à travers les murs en régime
   stationnaire en fonction de la profondeur fractale (0, 1, 2, 3).
