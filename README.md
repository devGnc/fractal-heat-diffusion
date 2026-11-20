# Diffusion thermique et géométries fractales

Projet de recherche personnel de Devrim Genc. Ce dépôt reprend la simulation
de l'ancien dépôt GitLab (`tests/diffusion`).

## Problématique

> Dans un monde qui se réchauffe, à aire égale, de combien une géométrie fractale
> augmente-t-elle la chaleur évacuée passivement, et à quelles conditions ?

Une enveloppe d'habitation a intérêt à être compacte : plus de surface, c'est plus de
pertes en hiver et plus d'apports en canicule. La fractale est donc étudiée comme
**échangeur** (façade, mur, dissipateur) chargé d'évacuer de la chaleur produite à
l'intérieur vers un milieu plus froid (air nocturne, sol).

On compare une pièce dont le contour est un **flocon de Koch** à une **pièce carrée
de même aire**, chauffées par une source centrale à 100 °C, avec l'extérieur maintenu
à 20 °C. À aire égale, la fractale a un périmètre bien plus long : elle échange donc
davantage de chaleur avec l'extérieur. La grandeur mesurée est la puissance
`P_source` évacuée en régime établi ; divisée par (T_source − T_ext), c'est la
conductance thermique de la pièce.

Mots-clés : fractales, suites, récursivité, programmation en C, équation de la
chaleur, différences finies, complexité.

## Compilation et exécution

```bash
brew install sdl2 pkg-config      # macOS
make                              # compile ./thermalfractalsim
make run                          # fenêtre interactive
make data                         # carré + Koch 0 à 3 sans fenêtre (≈ 15 s)
make plot                         # courbes temporelles + étude en profondeur
make test                         # validation : anneau comparé à la solution exacte (≈ 2 min)
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
boucle à la précision machine en régime établi (températures en `double`).

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
tests/test_anneau.c     validation sur un anneau (solution exacte), plusieurs résolutions
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
| Carré | 13 689 | 468 | 45,48 °C | 66,42 | 0,00 % |
| Koch 0 (triangle) | 13 795 | 664 | 42,90 °C | 69,46 | 0,00 % |
| Koch 1 | 13 757 | 764 | 42,72 °C | 69,59 | 0,00 % |
| Koch 2 | 13 771 | 960 | 41,47 °C | 71,27 | 0,00 % |
| Koch 3 | 13 645 | 1 188 | 41,13 °C | 71,96 | 0,00 % |

Valeurs calculées en `double`. Les CSV de `results/v2/` datent de la version en
`float` : P_source y est plus grande de 0,02 % et le bilan ne boucle qu'à −0,05 %.
Les pertes de Newton (`P_cool`) représentent environ 10 % de la puissance.

**Lecture.** De la profondeur 0 à 3, la surface d'échange augmente de 79 %, mais la
puissance seulement de 3,6 %. Avec ces paramètres, la résistance thermique est
dominée par la **conduction dans l'air immobile** entre la source et les murs, pas
par les murs : la surface gagnée au bout des pointes de la fractale, loin de la
source, sert peu.

Test avec un air 10 fois plus diffusif (α_air = 2, Δt = 0,1, ce qui imite un air
brassé par convection) : la puissance passe de 394 (Koch 0) à 471 (Koch 3), soit
**+20 %**, et Koch 3 dépasse le carré (360) de **+31 %**. La conclusion dépend donc
fortement du rapport α_air / α_mur, qu'il faut justifier physiquement.

## Validation : test de l'anneau

`make test` simule un anneau : source de rayon R1 = 8·s, air jusqu'au mur au rayon
32·s, extérieur au-delà, et compare la puissance à la solution exacte de la
conduction radiale, P = 2π ΔT / (ln(Ra/R1)/α_air + ln(Rb/Ra)/α_mur). Le code de la
simulation est réutilisé tel quel. Le cercle contient toutes les orientations de mur.
`R_mur` est la résistance du mur déduite de la simulation, rapportée à la valeur
exacte (1 = mur correctement représenté).

| Cas | Erreur sur P, s = 1 | s = 2 | s = 4 | R_mur (s = 4) |
|---|---|---|---|---|
| Milieu homogène (α_mur = α_air) | −3,6 % | −1,8 % | −0,8 % | — |
| Mur fin d'un pixel isolant (α_mur = α_air/100) | +12,3 % | +8,7 % | +5,9 % | 0,845 |
| Mur épais de 4·s pixels (α_mur = α_air/10) | −0,3 % | −0,1 % | −0,2 % | 1,004 |

**Lecture.**
- Le schéma converge à l'ordre 1 dans un milieu homogène : l'erreur est divisée par 2
  quand la résolution double. Elle vient de la position des bords en escalier
  (source et extérieur), connue à un pixel près.
- Un **mur d'un pixel** tracé comme les flocons (Bresenham) a une résistance trop
  faible de **15 à 16 %**, **quelle que soit la résolution** (0,843 ; 0,850 ; 0,845) :
  c'est le biais de grille du point 9, mesuré. L'erreur sur P baisse avec s
  seulement parce que le mur pèse de moins en moins dans la résistance totale.
- Un **mur épais** est bien représenté : sa résistance est juste à 0,7 % près dès
  s = 1, et à 0,4 % près à s = 4.
- Le cas « mur fin » avec les paramètres de la v2 (α_mur = α_air/10, non affiché)
  donne une erreur sur P très faible, mais par compensation : le biais du mur
  (−15 %) et l'erreur de l'air (+0,8 %, sur une résistance 18 fois plus grande)
  s'annulent presque. Une petite erreur sur un seul cas ne valide pas un modèle.
- En `float`, le régime établi se figeait avant d'être atteint (incréments plus
  petits que la précision) et le bilan s'écartait de 0,01 % à 0,16 %, d'autant plus
  que la grille était fine. Les températures sont désormais en `double`.

## Limites connues et pistes

1. ~~Objectif à préciser~~ : la fractale est étudiée comme échangeur qui évacue la
   chaleur (voir la problématique).
2. ~~Iso-aire approximative~~ : corrigé (écart de 0,6 %, au lieu de 4,6 %).
3. ~~Murs inclus dans la moyenne~~ : corrigé, la moyenne porte sur l'air hors source.
4. ~~Schéma non conservatif~~ : corrigé ; le bilan d'énergie boucle exactement
   depuis le passage en `double`.
5. ~~Régime non établi~~ : corrigé par le mode `--headless`.
6. ~~CFL très large~~ : Δt multiplié par 10.
7. **Unités de grille** : adimensionner ou fixer des Δx et Δt physiques.
8. **Paramètres physiques à justifier** :
   - le code met la diffusivité α là où il faudrait la conductivité λ
     (ρc·∂T/∂t = ∇·(λ∇T)), avec ρc = 1 partout. En régime établi seul λ compte :
     le béton conduit environ 70 fois mieux que l'air (λ), alors qu'il diffuse
     environ 25 fois moins (α). Le modèle a donc un mur 10 fois moins conducteur
     que l'air, l'inverse de la réalité ;
   - l'air d'une pièce n'est pas immobile : la convection le brasse, et la
     résistance se concentre près des parois (R_si ≈ 0,13 m²·K/W) et dans le mur.
     L'air immobile (α_air = 0,2) et des murs d'un pixel rendent les murs peu
     limitants ;
   - les pertes de Newton (`COOLING_FACTOR`, environ 10 % de la puissance)
     s'appliquent aussi aux murs ; à justifier (sol et plafond ?) ou à retirer. Pistes : air brassé (α_air effectif
   élevé ou température d'air homogène), murs plus épais, et tracer P_source en
   fonction de la profondeur pour plusieurs rapports α_air / α_mur.
9. **Biais de grille** : avec le stencil à 5 points, un mur oblique en escalier
   échange sur sa longueur |Δx|+|Δy| (flocon de profondeur 2 : 960 faces air/mur
   pour un périmètre géométrique d'environ 792, soit +21 %). Ce biais ne diminue pas en
   raffinant la grille : le test de l'anneau mesure une résistance de mur fin trop
   faible de 15 à 16 % à toutes les résolutions, alors qu'un mur épais de plusieurs
   pixels est juste à 0,4 % près. Comparer les profondeurs de Koch entre elles l'annule
   (mêmes orientations 0°, 60°, 120°) ; la comparaison avec le carré le subit.
