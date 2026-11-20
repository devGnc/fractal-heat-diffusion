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

.PHONY: all clean run plot data test

all: $(TARGET)

$(BUILDDIR):
	mkdir -p $(BUILDDIR)

$(TARGET): $(OBJ)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

$(BUILDDIR)/%.o: $(SRCDIR)/%.c $(HDR) | $(BUILDDIR)
	$(CC) $(CFLAGS) -c $< -o $@

run: $(TARGET)
	./$(TARGET)

# Simulations complètes sans fenêtre (régime établi), CSV à la racine :
# carré iso-aire et flocons de Koch de profondeur 0 à 3 (même aire)
GEOMETRIES = carre koch0 koch1 koch2 koch3
data: $(TARGET)
	@for g in $(GEOMETRIES); do ./$(TARGET) --headless $$g $(ITER) || exit 1; done

plot:
	python3 scripts/analyse.py
	python3 scripts/profondeur.py

# Test de validation : anneau circulaire comparé à la solution exacte, à
# plusieurs résolutions s (source de rayon 8·s, mur au rayon 32·s).
SCALES = 1 2 4
TEST_SRC = tests/test_anneau.c $(SRCDIR)/physics.c $(SRCDIR)/geometry.c
test: | $(BUILDDIR)
	@for s in $(SCALES); do \
		$(CC) $(CFLAGS) -DGRID_SIZE=$$((72 * s + 9)) -DSOURCE_RADIUS=$$((8 * s)) \
			-DCOOLING_FACTOR=0 -o $(BUILDDIR)/test_anneau_s$$s $(TEST_SRC) -lm || exit 1; \
	done
	@printf "%-13s %s %4s %8s %8s %11s %11s %8s %8s %9s %8s\n" \
		cas s grille Ra Rb P_simul P_exact "err_%" "bilan_%" R_mur Iter
	@for s in $(SCALES); do ./$(BUILDDIR)/test_anneau_s$$s || exit 1; done

# Ne supprime que les fichiers de compilation : les CSV de résultats sont conservés.
clean:
	rm -rf $(BUILDDIR) $(TARGET)
