Dynamic Genome Warping (DGW)
===============

# Installation
`numpy`, `scipy`, `pandas==0.10.1`, `mlpy`, `pysam`, `fastcluster` and `matplotlib` packages.


## MLPY
DGW depends on modified version of MLPY, available from https://github.com/sauliusl/mlpy
You can install it by
```
pip install cython
pip install -e git://github.com/sauliusl/mlpy.git#egg=mlpy
```
Cython is required for compilation of mlpy from source.
Note that mlpy depends on GSL (with header-files). Please look for instructions on how to install it on your platform if the installation fails.

## Platform-specific instructions

### Linux
You should be able to directly install DGW dependencies by using `pip` on linux:
```
pip install --uprgrade numpy scipy pandas mlpy pysam fastcluster matplotlib
```

Finally install modified MLPY:
```
pip install -e git://github.com/sauliusl/mlpy.git#egg=mlpy
```

### Mac OS X & Windows
It is probably easies to install `numpy`, `scipy`, `pandas` and `matplotlib` on Windows via Enthought python distribution. This distribution is free for academic use. See 
http://www.enthought.com/products/epd.php for instructions how to download and install it.

Unfortunately EPD 7.3 ships with an outdated version of `pandas` package that is not compatible with DGW, therefore you would have to upgrade it using pip in order to run it.
```
pip install --upgrade pandas
```

Packages that are not in EPD also need to be installed using pip
```
pip install pysam fastcluster
```

if `fastcluster` install fails install it from source from http://github.com/sauliusl/fastcluster

Don't forget MLPY:
```
pip install -e git://github.com/sauliusl/mlpy.git#egg=mlpy
```

