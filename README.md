# EM-BDD Algorithm Implementation

This repository contains the complete source code for the EM-BDD (Expectation-Maximization learning with Binary Decision Diagrams) algorithm developed for learning Hidden Markov Models (HMMs). 
The EM-BDD algorithm leverages Binary Decision Diagrams to optimize the learning process, offering a alternative to traditional methods such as the Baum-Welch algorithm.
Highly inspired by the work of Ishihata and colleagues in the paper [An EM algorithm on BDDs with order encoding for logic-based probabilistic models](https://proceedings.mlr.press/v13/ishihata10a).


For more details on the EM-BDD algorithm developed here, see the official publication: [The EM-BDD Algorithm for Learning Hidden Markov Models](https://link.springer.com/chapter/10.1007/978-3-031-75107-3_8).

## Authors
- [Eva Ósk Gunnarsdóttir](https://github.com/efaosk)
- [Anna Ingólfsdóttir](https://scholar.google.com/citations?user=B0YC1p8AAAAJ&hl=en&oi=ao)

## Overview
The EM-BDD algorithm aims to enhance the efficiency of HMM learning by using Binary Decision Diagrams (BDDs) to represent the state space. This implementation provides both EM-BDD and traditional Baum-Welch algorithms to enable comparisons and assess improvements in training Hidden Markov Models.

## Features

- **EM-BDD Algorithm**: A novel approach for HMM training that utilizes BDDs for efficient probability computation and state representation.
- **Baum-Welch (BW) Algorithm**: The standard algorithm for HMM training, included for comparative analysis.

## Getting Started

### Prerequisites

- **C Compiler**: Ensure you have GCC or another compatible C compiler installed.
- **CUDD Library**: This project requires the CUDD library for BDD manipulation. Install it following the instructions on the [CUDD website](https://davidkebo.com/cudd/)

To set up the CUDD library paths, add the include and lib directories to CPATH and LIBRARY_PATH, respectively:
```bash
export CPATH=/path/to/cudd/includes
export LIBRARY_PATH=/path/to/cudd/libs
```

### Installation

1. Clone the repository to your local machine:

    ```bash
    git clone https://github.com/EfaOsk/EM-BDD-Algorithm.git
    cd EM-BDD
    ```

2. Compile the Project using the provided Makefile:
    ```bash
    make
    ```

3. Run the Program:
    ```bash
    ./main <dataset_file> <num_states> <epsilon>
    ```

    - `dataset_file`: Path to the dataset file containing observation sequences.
    - `num_states`: Number of states in the HMM model.
    - `epsilon`: Convergence threshold for the EM algorithm.

    Example:
    ```bash
    ./main example_dataset.txt 3 0.01
    ```

### Dataset Format
The dataset file should be formatted as follows:

- First line: Number of possible observations ($M$).
- Second line: Length of the observation sequences ($T$).
- Third line: Number of sequences in the dataset ($K$).
- Subsequent lines: Each line contains a sequence of observations, separated by spaces.

See file `example_dataset.txt`.

### Project Structure
- `main.c`: The main file that loads the dataset, runs the EM-BDD and Baum-Welch algorithms, and logs results.
- `helpers.c` / `helpers.h`: Utility functions for data handling, matrix operations, and BDD utilities.
- `BDD_build.c`: Functions for creating and manipulating Binary Decision Diagrams.
- `EMBDD.c`: Implementation of the EM-BDD algorithm.
- `HMM_algorithms.c`: Standard HMM operations.
- `HMM_management.c`: HMM model management and initialization.
- `exampleHMM.c`: Sample HMM initializations for testing and experimentation.