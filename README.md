# Orthogonal Decompositions and the Matchgate Group
This repository contains Python functions designed to implement the Hoffman algorithm for decomposing elements of SO(N)[1,2], and the application of this algorithm to decomposing any element of the matchgate group into nearest-neighbor matchgates[3]. The Hoffman algorithm may be of independent interest. The decomposition of matchgate group elements is a key ingredient in the matchgate randomized benchmarking introduced in [4], and is also the most efficient method to generate elements of the matchgate group distributed according to the Haar measure. We refer to [2] for details on the Hoffman algorithm, and [4] for details on the matchgate group, matchgate decompositions, and matchgate benchmarking. The code has two sections.

## Hoffman Decomposition
First, it implements the Hoffman decomposition of matrices in SO(N) into N(N-1)/2 rotations involving two basis elements each, first described in Hoffman, Raffenetti, and Ruedenberg, J. Math. Phys. 13, 528 (1972). This is done in the file SpecialOrthogonalDecomposition.py.

## Matchgate Circuits
Second, it implements the decomposition of an arbitrary element of SO(2N) into unitary matchgates acting on nearest neighbors. This decomposition was first described in Jozsa and Miyake, Proc. Math. Phys. Eng. Sci. 464, 3089 (2008), although we use a slightly modified decomposition that requires fewer matchgates. This is done in the file MatchgateDecomposition.py, which requires the functions in SpecialOrthogonalDecomposition.py.

# References
[1] R. C. Raffenetti and K. Ruedenberg, Parametrization of an orthogonal matrix in terms of generalized Eulerian angles, Int. J. Quantum Chem. 4, 625 (1969).
[2] D. K. Hoffman, R. C. Raffenetti, and K. Ruedenberg, Generalization of Euler angles to N-dimensional orthogonal matrices, J. Math. Phys. 13, 528 (1972).
[3] R. Jozsa and A. Miyake, Matchgates and classical simulation of quantum circuits, Proc. Math. Phys. Eng. Sci. 464, 3089 (2008).
[4] Forthcoming
