#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

#include "HMM.h"
#include "learn.h"

int main(int argc, char *argv[]) {
    // Check for the correct number of arguments
    if (argc != 4) {
        fprintf(stderr, "Usage: %s N M T\n where:\n\tN = number of states\n\tM = alphabet size\n\tT = length of observation sequence\n", argv[0]);
        return 1; // Return an error code
    }

    // Parse command-line arguments and validate them
    int N, M, T;
    if (sscanf(argv[1], "%d", &N) != 1 || sscanf(argv[2], "%d", &M) != 1 || sscanf(argv[3], "%d", &T) != 1) {
        fprintf(stderr, "Invalid arguments. N, M, and T must be integers.\n");
        return 1; // Return an error code
    }

    // The observation O
    int O[T];
    O[0] = 0;
    O[1] = 1;

    // Call the learn function
    struct HMM L = learn(N, M, T, O);

    // Clean up and return
    return 0;
}
