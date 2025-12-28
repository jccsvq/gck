# -----------------------------------------------------------------------------
# Makefile to compile the geostatistical programs in ./bin
# -----------------------------------------------------------------------------

# --- Configuration Variables ---

# Directory for the executables
BIN_DIR = bin

# Compiler and Compilation Flags
CC = cc
CFLAGS = -Wall -Wextra -O2 

# Common functions library and linking options
LIBRARY_SOURCE = libGCK2.c
LDLIBS = -lm

# Program Names (Executables)
PROGRAM_NAMES = gik2 geko crosval2 faik2

# Main source files containing main()
SOURCES = $(PROGRAM_NAMES:=.c)

# Full paths of the executables (bin/gik2, bin/geko, etc.)
PROGRAMS = $(addprefix $(BIN_DIR)/, $(PROGRAM_NAMES))

# Object files (created in the current directory)
OBJECTS = $(SOURCES:.c=.o)
LIBRARY_OBJECT = $(LIBRARY_SOURCE:.c=.o)

# -----------------------------------------------------------------------------
# 1. Main Rule ('all' Target)
# -----------------------------------------------------------------------------
# The main goal is to build all programs.
# It depends on the creation of the bin directory and all the executables.
.PHONY: all
all: $(BIN_DIR) $(PROGRAMS)

# -----------------------------------------------------------------------------
# 2. Target to create the 'bin' directory (if it does not exist)
# -----------------------------------------------------------------------------
# Ensures that '...'
# The 'mkdir -p' operator creates the directory if it does not exist.
$(BIN_DIR):
	@mkdir -p $(BIN_DIR)

# -----------------------------------------------------------------------------
# 3. Shared Library Compilation
# -----------------------------------------------------------------------------
# Rule to compile libGCK2.c to libGCK2.o (Library Object)
$(LIBRARY_OBJECT): $(LIBRARY_SOURCE)
	@echo "Compiling the library $(LIBRARY_SOURCE) to $(LIBRARY_OBJECT)..."
	$(CC) $(CFLAGS) -c $< -o $@

# -----------------------------------------------------------------------------
# 4. Compilation and Linking Rules
# -----------------------------------------------------------------------------

# Rule to compile each source file to its object (.o)
%.o: %.c
	@echo "Compiling $^ to $@"
	$(CC) $(CFLAGS) -c $< -o $@

# Linking Rule: Moves the executables to $(BIN_DIR)
# Depends on the main object (.o) and the library object (libGCK2.o)
# The vertical bar syntax (|) ensures that $(BIN_DIR) is executed first.
$(BIN_DIR)/%: %.o $(LIBRARY_OBJECT) | $(BIN_DIR)
	@echo "Linking $@..."
	$(CC) $(CFLAGS) $< $(LIBRARY_OBJECT) $(LDLIBS) -o $@

# -----------------------------------------------------------------------------
# 5. Clean Target
# -----------------------------------------------------------------------------
.PHONY: clean
clean:
	@echo "Cleaning object files and executables..."
	# Delete temporary objects
	rm -f $(OBJECTS) $(LIBRARY_OBJECT)
	# Delete executables in bin
	rm -f $(PROGRAMS)
	# Optional: to remove the bin directory if it is empty
	@-rmdir $(BIN_DIR) 2>/dev/null || true
