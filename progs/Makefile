# -----------------------------------------------------------------------------
# Makefile to compile C programs into the 'bin/' directory
# -----------------------------------------------------------------------------

# The target directory for the executables
BIN_DIR = bin

# The C compiler we will use
CC = gcc

# Compiler flags.
CFLAGS = -Wall -Wextra -O2

# Libraries to link. We need the math library.
LDLIBS = -lm

# List of all executable names (without path)
PROGRAM_NAMES = gik2 geko crosval2 faik2

# List of all executables with the BIN_DIR/ path
# Use the 'addprefix' Make function to add "bin/" to each name.
PROGRAMS = $(addprefix $(BIN_DIR)/, $(PROGRAM_NAMES))

# Source files (in the current directory)
SOURCES = $(PROGRAM_NAMES:=.c)

# -----------------------------------------------------------------------------
# Default target: Compiles all programs
# -----------------------------------------------------------------------------
.PHONY: all
all: $(BIN_DIR) $(PROGRAMS)

# -----------------------------------------------------------------------------
# Target to create the 'bin' directory
# -----------------------------------------------------------------------------
# This target depends on the $(BIN_DIR) directory existing.
# The 'mkdir -p' operator creates the directory if it does not exist.
$(BIN_DIR):
	@mkdir -p $(BIN_DIR)

# -----------------------------------------------------------------------------
# Compilation rules
# -----------------------------------------------------------------------------

# Generic rule to compile a .c file to an executable in the bin directory.
# Pre-requisite: $(BIN_DIR) to ensure the directory exists before compiling.
# $@: the target name (e.g., bin/gik2)
# $*: the part of the name that matches the '%' pattern (e.g., gik2)
# We use the prefix expansion of the target to get the source file
$(BIN_DIR)/%: %.c | $(BIN_DIR)
	@echo "Compiling $< to $@"
	$(CC) $(CFLAGS) $< -o $@ $(LDLIBS)

# -----------------------------------------------------------------------------
# Explicit rules so the 'all' target knows what to look for (for clarity)
# -----------------------------------------------------------------------------

$(BIN_DIR)/gik2: gik2.c
$(BIN_DIR)/geko: geko.c
$(BIN_DIR)/crosval2: crosval2.c
$(BIN_DIR)/faik2: faik2.c

# -----------------------------------------------------------------------------
# clean target: Delete executables
# -----------------------------------------------------------------------------
.PHONY: clean
clean:
	@echo "Cleaning executable files in $(BIN_DIR)..."
	rm -f $(PROGRAMS)
	# Optional: to remove the bin directory if it is empty
	@-rmdir $(BIN_DIR) 2>/dev/null || true
