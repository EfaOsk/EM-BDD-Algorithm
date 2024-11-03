# Compiler and flags
CC = gcc
CFLAGS = -I$(CPATH) -Icudd -Iutil -lm

# Path to CUDD library
CUDD_LIB ?= $(LIBRARY_PATH)/libcudd.a

# Output executable
TARGET = main

# Source files
SRCS = helpers.c BDD_build.c EMBDD.c HMM_utils.c HMM_algorithms.c HMM_management.c exampleHMM.c main.c

# Object files (each source file will have a corresponding .o file)
OBJS = $(SRCS:.c=.o)

# Default target: build the executable
all: $(TARGET)

# Link all object files to create the final executable
$(TARGET): $(OBJS)
	$(CC) -o $(TARGET) $(OBJS) $(CUDD_LIB) $(CFLAGS)

# Compile .c files to .o files
%.o: %.c
	$(CC) -c $< -o $@ $(CFLAGS)

# Run the executable with default arguments
run: $(TARGET)
	./$(TARGET) example_dataset.txt 3 0.01

# Clean up build artifacts
clean:
	rm -f $(TARGET) $(OBJS)

# Rebuild everything from scratch
rebuild: clean all
