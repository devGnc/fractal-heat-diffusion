# Nom de l'exécutable final
TARGET = thermalfractalsim

# Répertoires
SRCDIR = src
BUILDDIR = build

# Compilateur et options
CC = gcc
CFLAGS = -Wall -Wextra -O2 -std=c11 `pkg-config --cflags sdl2` -I$(SRCDIR)
LDFLAGS = `pkg-config --libs sdl2` -lm

# Fichiers sources (dans src/)
SRC = $(wildcard $(SRCDIR)/*.c)
HDR = $(wildcard $(SRCDIR)/*.h)

# Fichiers objets (dans build/)
OBJ = $(patsubst $(SRCDIR)/%.c,$(BUILDDIR)/%.o,$(SRC))

.PHONY: all clean run plot data regime loi test rapport animation reencoder figures

all: $(TARGET)

$(BUILDDIR):
	mkdir -p $(BUILDDIR)

$(TARGET): $(OBJ)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

$(BUILDDIR)/%.o: $(SRCDIR)/%.c $(HDR) | $(BUILDDIR)
	$(CC) $(CFLAGS) -c $< -o $@

run: $(TARGET)
	./$(TARGET)

GEOMETRIES = carre koch0 koch1 koch2 koch3 koch4

# Courbes dans le temps, sans fenêtre, jusqu'au régime établi (grille de la
# fenêtre, Δx = 3,5 cm ; ≈ 2 min) : res_<geometrie>.csv à la racine.
data: $(TARGET)
	@for g in $(GEOMETRIES); do ./$(TARGET) --headless $$g || exit 1; done

# Régime établi précis, résolu directement par gradient conjugué sur une
# grille fine (Δx = 1,25 cm ; ≈ 1 min) : regime_<geometrie>.csv à la racine.
FINE = $(BUILDDIR)/thermalfractalsim_fin
$(FINE): $(SRC) $(HDR) | $(BUILDDIR)
	$(CC) $(CFLAGS) -DGRID_SIZE=1200 "-DDX=0.0125" -o $@ $(SRC) $(LDFLAGS)

regime: $(FINE)
	@for g in $(GEOMETRIES); do ./$(FINE) --regime $$g || exit 1; done

# Règle empirique de la profondeur optimale : régime établi sur la grille fine pour 21
# épaisseurs de mur (5 à 45 cm), 36 m², flocons, carrés crénelés et îles de
# Minkowski (273 calculs, tous présents dans results/v4/loi/loi.csv : la cible
# ne recalcule que ce qui manque).
loi: $(FINE)
	python3 scripts/balayage.py

# Animations de la fenêtre, filmée sans l'ouvrir (une semaine simulée, ≈ 1 min) :
# GIF léger pour le README (640 px, 10 images/s, 64 couleurs, ≈ 1,5 Mo), MP4 et
# WebM pour le web (840 px, quelques centaines de Ko).
ANIM = docs/captures/diffusion
ANIM_FILTRES = [0]split=3[g][m][w];\
	[g]fps=10,scale=640:-1:flags=lanczos,split[a][b];[a]palettegen=max_colors=64[p];\
	[b][p]paletteuse=dither=bayer:bayer_scale=5:diff_mode=rectangle[gif];\
	[m]scale=840:-2:flags=lanczos,format=yuv420p[mp4];[w]scale=840:-2:flags=lanczos[webm]
ANIM_SORTIES = -map "[gif]" $(ANIM).gif \
	-map "[mp4]" -c:v libx264 -preset veryslow -crf 30 -movflags +faststart -an $(ANIM).mp4 \
	-map "[webm]" -c:v libvpx-vp9 -crf 40 -b:v 0 -row-mt 1 -an $(ANIM).webm
animation: $(TARGET)
	./$(TARGET) --film 3 200 | ffmpeg -loglevel error -y -f rawvideo -pixel_format bgra \
		-video_size 1400x780 -framerate 20 -i - \
		-filter_complex "[0]tpad=stop_mode=clone:stop_duration=2[s];$(subst [0],[s],$(ANIM_FILTRES))" \
		$(ANIM_SORTIES)
	@ls -lh $(ANIM).*

# Réencode les animations à partir d'un GIF déjà produit, sans simulation :
# make reencoder SOURCE=chemin/vers/source.gif
reencoder:
	ffmpeg -loglevel error -y -i $(SOURCE) -filter_complex "$(ANIM_FILTRES)" $(ANIM_SORTIES)
	@ls -lh $(ANIM).*

plot:
	python3 scripts/analyse.py
	python3 scripts/profondeur.py
	python3 scripts/loi.py

# Test de validation : pièce circulaire comparée à la solution exacte, à
# plusieurs résolutions s (Δx = 5 cm / s).
SCALES = 1 2 4
TEST_SRC = tests/test_anneau.c $(SRCDIR)/physics.c $(SRCDIR)/geometry.c $(SRCDIR)/io.c
test: | $(BUILDDIR)
	@for s in $(SCALES); do \
		$(CC) $(CFLAGS) -DGRID_SIZE=$$((80 * s + 1)) "-DDX=(0.05/$$s)" \
			-o $(BUILDDIR)/test_anneau_s$$s $(TEST_SRC) -lm || exit 1; \
	done
	@printf "%-12s %s %4s %9s %9s %9s %9s %8s %10s %8s\n" \
		cas s grille T_air T_exacte G_simul G_exact "err_G_%" "bilan_%" Iter
	@for s in $(SCALES); do ./$(BUILDDIR)/test_anneau_s$$s || exit 1; done

# Figures du rapport et du README (docs/rapport/figures/), à partir de results/v4
figures:
	python3 docs/rapport/figures/generer_figures.py

# Rapport LaTeX (docs/rapport/rapport.pdf)
rapport:
	cd docs/rapport && latexmk -pdf -interaction=nonstopmode -halt-on-error rapport.tex

# Ne supprime que les fichiers de compilation : les CSV de résultats sont conservés.
clean:
	rm -rf $(BUILDDIR) $(TARGET)
