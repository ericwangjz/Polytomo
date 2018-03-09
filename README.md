# Polytomo

Overview
--------
This toolbox aims to compute confidence region in quantum state tomography.

It consists of three modules:
1. polyconfiregion
2. sampling
3. polytomorun

The `polyconfiregion` module consists of tools that contruct the polytope confidence region
The `sampling` module draws confidence intervals from the polytope via sampling
The `polytomorun` is the interface that reads all the input information, such as measurements POVM and respective counts, and outputs the figures of merit.

Prerequisites
------

1. the Python `tomographer` package by P. Faist  available at https://tomographer.github.io/tomographer/
2. the `polytope 0.2.1` package available at https://pypi.python.org/pypi/polytope/0.2.1
3. `qutip` package available at http://qutip.org/docs/latest/index.html



