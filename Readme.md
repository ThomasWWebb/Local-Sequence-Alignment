# Bioinformatics - Local sequence alignment

Implementations of three algorithms for local sequence alignment. Two dynamic approaches; `dynprog` that runs in quadratic time and space, and `dynproglin` that runs in linear space. The third is an impelmentation of FASTA `heuralign`, a heuristic algorithm that runs in sub-quadratic time

Each function takes the same 4 inputs: a string string of all (unique) letters of length p, a symmetric (p + 1) x (p + 1) substitution matrix (represented by a list of lists) and the two sequences to align.
