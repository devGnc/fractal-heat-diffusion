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

.PHONY: all clean run plot

all: $(TARGET)

$(BUILDDIR):
	mkdir -p $(BUILDDIR)

$(TARGET): $(OBJ)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

$(BUILDDIR)/%.o: $(SRCDIR)/%.c $(HDR) | $(BUILDDIR)
	$(CC) $(CFLAGS) -c $< -o $@

run: $(TARGET)
	./$(TARGET)

plot:
	python3 scripts/analyse.py

# Ne supprime que les fichiers de compilation : les CSV de résultats sont conservés.
clean:
	rm -rf $(BUILDDIR) $(TARGET)
