pycute - Pythonic wrapper for Correlation Utilities and Two-point Estimates
===========================================================================

Introduction
------------

**pycute** calculates 2-pt, 3-pt (following https://arxiv.org/abs/1506.02040) and 4-pt functions required for integral constraint corrections (https://arxiv.org/abs/1904.08851v3).
This package is inspired by CUTE https://github.com/damonge/CUTE (https://arxiv.org/abs/1210.1833).
The C code has been mostly refurbished into a more flexible (though probably slower) form. It uses OpenMP for parallelization if available.
A basic Python wrapper (based on ctypes) is provided in pycute.py.
Different line-of-sight definitions are implemented (endpoint, midpoint).
All tests are in tests/.

Installation
------------

To compile the C code:\
$ make clean\
$ make

Requirements
------------

- scipy
- numpy
