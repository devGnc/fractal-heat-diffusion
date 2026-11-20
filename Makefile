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

# Simulations complètes sans fenêtre jusqu'au régime établi, CSV à la racine :
# carré et flocons de Koch de profondeur 0 à 4 (même aire). Grille fine
# (Δx = 1,25 cm) : à 5 cm, les détails de profondeur 3 et 4 sont mal résolus.
FINE = $(BUILDDIR)/thermalfractalsim_fin
$(FINE): $(SRC) $(HDR) | $(BUILDDIR)
	$(CC) $(CFLAGS) -DGRID_SIZE=1200 "-DDX=0.0125" -o $@ $(SRC) $(LDFLAGS)

GEOMETRIES = carre koch0 koch1 koch2 koch3 koch4
data: $(FINE)
	@for g in $(GEOMETRIES); do ./$(FINE) --headless $$g || exit 1; done

plot:
	python3 scripts/analyse.py
	python3 scripts/profondeur.py

# Test de validation : pièce circulaire comparée à la solution exacte, à
# plusieurs résolutions s (Δx = 5 cm / s ; Δt réduit en 1/s² pour la stabilité).
SCALES = 1 2 4
TEST_SRC = tests/test_anneau.c $(SRCDIR)/physics.c $(SRCDIR)/geometry.c
test: | $(BUILDDIR)
	@for s in $(SCALES); do \
		$(CC) $(CFLAGS) -DGRID_SIZE=$$((80 * s + 1)) "-DDX=(0.05/$$s)" \
			-DDT=$$((200 / (s * s))).0 -o $(BUILDDIR)/test_anneau_s$$s $(TEST_SRC) -lm || exit 1; \
	done
	@printf "%-12s %s %4s %10s %10s %8s %9s %8s %8s\n" \
		cas s grille P_simul P_exact "err_%" "bilan_%" R_surf Iter
	@for s in $(SCALES); do ./$(BUILDDIR)/test_anneau_s$$s || exit 1; done

# Ne supprime que les fichiers de compilation : les CSV de résultats sont conservés.
clean:
	rm -rf $(BUILDDIR) $(TARGET)
