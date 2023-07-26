#!/bin/bash

# Define the values of N and M
N_values=(2 5 10 15)
M_values=(2 5 10 15)

# Loop over the values of N and M
for N in "${N_values[@]}"
do
    for M in "${M_values[@]}"
    do
        # Run the ./test program and save the output to results.txt
        ./test $N $M >> results.txt
    done
done
