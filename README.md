# TIPE : diffusion thermique et géométries fractales

TIPE MPI de Devrim Genc. Ce dépôt reprend la simulation de
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
make data                         # carré + Koch 0 à 3 sans fenêtre (≈ 15 s)
make plot                         # courbes temporelles + étude en profondeur
```

Mode sans fenêtre, pour obtenir des données exploitables (enregistrement dès
l'itération 0, une ligne toutes les 400 itérations) :

```bash
./thermalfractalsim --headless carre                # 60 000 itérations → res_carre.csv
./thermalfractalsim --headless koch3                # flocon de profondeur 3 → res_koch3.csv
./thermalfractalsim --headless koch2 100000 f.csv   # nombre d'itérations et fichier choisis
make data ITER=100000
```

Géométries : `carre`, `koch0` (triangle) à `koch4`, et `fractal` (= `koch2`).
Le résumé final affiche le **bilan d'énergie** `P_source = P_ext + P_cool`, qui
doit boucler à mieux que 0,1 % en régime établi.

| Touche | Action |
|---|---|
| Espace | Basculer carré ↔ fractale |
| 0 à 4 | Flocon de Koch de cette profondeur |
| S | Démarrer / arrêter l'enregistrement CSV (`res_<geometrie>.csv`) ; changer de géométrie ou R l'arrête aussi |
| R | Réinitialiser |
| Q | Quitter |

Format du CSV : une ligne de commentaire décrivant la géométrie
(`# geometrie=koch2 profondeur=2 aire=... faces_air_mur=...`, lue avec
`pandas.read_csv(..., comment="#")`), puis les colonnes
`Iter,Time,Mean,Std,Min,P_source,P_ext,P_cool` :

- `Mean`, `Std`, `Min` : température de l'**air intérieur**, murs et source exclus.
- `P_source` : puissance fournie par le radiateur pour rester à 100 °C, c'est-à-dire
  l'énergie réinjectée dans la source à chaque pas divisée par Δt (°C·px par unité
  de temps, capacité thermique 1 par pixel).
- `P_ext` : puissance sortant par les faces mur/extérieur ; `P_cool` : pertes de
  Newton. En régime établi, `P_source = P_ext + P_cool`.

Les CSV produits à la racine sont ignorés par git : il faut les déplacer dans
`results/<version>/` pour les garder.

## Organisation

```
src/
  thermal_diffusion.h   paramètres physiques et structure Sim
  geometry.c            Koch récursif iso-aire, Bresenham, flood fill, carré iso-aire
  physics.c             schéma conservatif, conductances des faces, CFL, statistiques, bilan
  render.c              carte de couleurs SDL
  main.c                mode headless, boucle SDL, clavier, enregistrement CSV
scripts/analyse.py      évolution temporelle de toutes les simulations d'un dossier
scripts/profondeur.py   régime établi en fonction de la profondeur de Koch
results/v0/             premiers résultats (ancien dépôt, régime non établi)
results/v1/             headless, stats sur l'air, schéma encore non conservatif
results/v2/             schéma conservatif, carré + Koch 0 à 3
docs/                   export GeoGebra
```

## Modèle actuel (v2)

- Grille 300×300, Δx = 1, **Δt = 0,5**, α_air = 0,20, α_mur = 0,02 ; murs d'un pixel
- **Schéma conservatif** (volumes finis) : chaque face a une conductance égale à la
  moyenne harmonique des α des deux cases (deux demi-cases en série). La face
  mur/extérieur vaut 2·α_mur : l'extérieur est un thermostat parfait à 20 °C.
- Stabilité : Δt·Σk ≤ 1 pour chaque case (ici 0,4)
- Source : disque de rayon 15 à 100 °C ; pertes de Newton faibles (`COOLING_FACTOR`)
- Flocons de profondeur 0 à 3 de **même aire théorique**
  (A_n = A_0 · (8/5 − 3/5 · (4/9)^n), le triangle de départ est mis à l'échelle) ;
  carré iso-aire à 0,6 % près du flocon de profondeur 2

`make data` (régime établi vers t ≈ 15 000, soit 30 000 itérations) :

| | Aire (px) | Faces air/mur | T moy. air | P_source | Bilan |
|---|---|---|---|---|---|
| Carré | 13 689 | 468 | 45,48 °C | 66,44 | −0,06 % |
| Koch 0 (triangle) | 13 795 | 664 | 42,90 °C | 69,48 | −0,05 % |
| Koch 1 | 13 757 | 764 | 42,72 °C | 69,60 | −0,05 % |
| Koch 2 | 13 771 | 960 | 41,47 °C | 71,28 | −0,05 % |
| Koch 3 | 13 645 | 1 188 | 41,13 °C | 71,97 | −0,05 % |

**Lecture.** De la profondeur 0 à 3, la surface d'échange augmente de 79 %, mais la
puissance seulement de 3,6 %. Avec ces paramètres, la résistance thermique est
dominée par la **conduction dans l'air immobile** entre la source et les murs, pas
par les murs : la surface gagnée au bout des pointes de la fractale, loin de la
source, sert peu.

Test avec un air 10 fois plus diffusif (α_air = 2, Δt = 0,1, ce qui imite un air
brassé par convection) : la puissance passe de 394 (Koch 0) à 471 (Koch 3), soit
**+20 %**, et Koch 3 dépasse le carré (360) de **+31 %**. La conclusion dépend donc
fortement du rapport α_air / α_mur, qu'il faut justifier physiquement.

## Limites connues et pistes

1. **Objectif à préciser** : la fractale isole moins bien mais dissipe mieux
   (radiateurs, dissipateurs, façades qui évacuent la chaleur).
2. ~~Iso-aire approximative~~ : corrigé (écart de 0,6 %, au lieu de 4,6 %).
3. ~~Murs inclus dans la moyenne~~ : corrigé, la moyenne porte sur l'air hors source.
4. ~~Schéma non conservatif~~ : corrigé ; le bilan d'énergie boucle à 0,05 %.
5. ~~Régime non établi~~ : corrigé par le mode `--headless`.
6. ~~CFL très large~~ : Δt multiplié par 10.
7. **Unités de grille** : adimensionner ou fixer des Δx et Δt physiques.
8. **Paramètres physiques à justifier** : l'air immobile (α_air = 0,2) et des murs
   d'un pixel rendent les murs peu limitants. Pistes : air brassé (α_air effectif
   élevé ou température d'air homogène), murs plus épais, et tracer P_source en
   fonction de la profondeur pour plusieurs rapports α_air / α_mur.
9. **Biais de grille** : avec le stencil à 5 points, un mur oblique en escalier
   échange sur sa longueur |Δx|+|Δy| (flocon de profondeur 2 : 960 faces air/mur
   pour un périmètre géométrique d'environ 792, soit +21 %). Ce biais ne diminue pas en
   raffinant la grille. Comparer les profondeurs de Koch entre elles l'annule
   (mêmes orientations 0°, 60°, 120°) ; la comparaison avec le carré le subit.
