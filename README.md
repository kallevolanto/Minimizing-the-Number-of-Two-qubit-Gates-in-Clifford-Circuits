### This repository contains the Matlab code used for calculating the number of two-qubit gates needed for the two algorithms discussed in the thesis and the code generating the uniformly random symplectic matrices.

### 'F = rand_symp_mat(n)' to generate a random symplectic matrix acting on n qubits

### '[N2, N1, NSwaps] = two_qubit_Tv_decomp(F)' to get the number of 2-qubit transvections, 1-qubit transvections and swap gates obtained using the 2-qubit transvections decomposition algorithm

### 'N2 = bruhat_gate_decomposition(F,3)' to get the number of 2-qubit gates obtained using the Bruhat decomposition algorithm
