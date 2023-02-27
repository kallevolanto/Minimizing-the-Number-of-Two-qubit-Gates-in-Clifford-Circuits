This repository contains the Matlab code used for calculating the number of two-qubit gates needed for the two algorithms discussed in the thesis and the code generating the uniformly random symplectic matrices.







#
### 'F = rand_symp_mat(n)' to generate a random symplectic matrix acting on n qubits

Can, Trung. "An algorithm to generate a unitary transformation from logarithmically many random bits." Research Independent Study (2017).








#
### '[N2, N1, NSwaps] = two_qubit_Tv_decomp(F)' to get the number of 2-qubit transvections, 1-qubit transvections and swap gates obtained using the 2-qubit transvections decomposition algorithm









#
### 'N2 = bruhat_gate_decomposition(F,3)' to get the number of 2-qubit gates obtained using the Bruhat decomposition algorithm

D. Maslov and M. Roetteler, “Shorter stabilizer circuits via bruhat decomposition and quantum circuit transformations,” IEEE Transactions on Information Theory, vol. 64, no. 7, pp. 4729–4738, 2018.

S. Aaronson and D. Gottesman, “Improved simulation of stabilizer circuits, ”Physical Review A, vol. 70, no. 5, p. 052328, 2004.

G. Strang, “Fast transforms: Banded matrices with banded inverses,” Proceedings of the National Academy of Sciences, vol. 107, no. 28, pp. 12413–12416, 2010.
