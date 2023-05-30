# Hoffman Decomposition and the Matchgate Group
This repository contains Python functions implementing the Hoffman algorithm for decomposing elements of SO(N)[1,2], and functions applying this algorithm to decompose any element of the matchgate group into nearest-neighbor matchgates[3]. Although we use the Hoffman algorithm here exclusively for decomposing matchgates, it may be of independent interest. The decomposition of matchgate group elements is a key ingredient in the matchgate randomized benchmarking introduced in [4], and is also an efficient method to generate elements of the matchgate group distributed according to the Haar measure on SO(2N). We refer to [2] for details on the Hoffman algorithm, and [4] for details on the matchgate group, matchgate decompositions, and matchgate benchmarking. The code has two sections.

## Hoffman Decomposition
The file HoffmanDecomposition.py contains functions that implement the iterative Hoffman algorithm to decompose any matrix in SO(N) into the product of N(N-1)/2 rotations, each of which only acts on two basis vectors[1,2]. The algorithm parameterizes every matrix T in SO(N) by N(N-1)/2 indepdent Euler angles, and uses these angles to generate the N(N-1)/2 rotations. We provide functions to generate the Euler angles corresponding to a given T, as well as to generate T from the Euler angles. We also provide a function that realizes the Hoffman decomposition, using the Euler angles to output a list of matrices, each of which act nontrivially on only 2 of the N dimensions, and which compose to T.

## Matchgate Circuits
The file MatchgateDecomposition.py realizes the bijection between SO(2N) and the matchgate group. Specifically, we provide a function that takes an arbitrary element T of SO(2N) and provides a list of less than 4N^3 nearest-neighbor matchgates that compose to implement T. This method was first introduced in [3], but we implement the simpler modification described in [4] for the SWAP operations. We also provide a function that takes a list of matchgates and outputs the element of SO(2N) corresponding to the matchgate circuit, fully realizing the bijection in both directions.

# References
[1] R. C. Raffenetti and K. Ruedenberg, Parametrization of an orthogonal matrix in terms of generalized Eulerian angles, Int. J. Quantum Chem. 4, 625 (1969).

[2] D. K. Hoffman, R. C. Raffenetti, and K. Ruedenberg, Generalization of Euler angles to N-dimensional orthogonal matrices, J. Math. Phys. 13, 528 (1972).

[3] R. Jozsa and A. Miyake, Matchgates and classical simulation of quantum circuits, Proc. Math. Phys. Eng. Sci. 464, 3089 (2008).

[4] J. Claes, E. Rieffel, and Z. Wang, Character randomized benchmarking for non-multiplicity-free groups with applications to subspace, leakage, and matchgate randomized benchmarking, PRX Quantum 2, 010351 (2021)
