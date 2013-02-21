Dynamic Genome Warping (DGW)
===============

Installation
============================
`numpy`, `scipy`, `pandas==0.10.1`, `mlpy`, `pysam`, `fastcluster` and `matplotlib` packages.

MLPY
======================
DGW depends on modified version of MLPY, available from https://github.com/sauliusl/mlpy
You can install it by
```
pip install -e git://github.com/sauliusl/mlpy.git#egg=mlpy
```

Suggested way: Enthought Distribution
=======================
`numpy`, `scipy`, `pandas` and `matplotlib` can be installed from Enthought distribution.


After it was installed, I recommend creating a `virtualenv` for it:
```
mkvirtualenv epd --python=/Library/Frameworks/EPD64.framework/Versions/7.3/bin/python --system-site-packages
```

note that this requires virtualenvwrapper: http://www.doughellmann.com/projects/virtualenvwrapper/

Unfortunately EPD 7.3 ships with an outdated version of `pandas` package that is not compatible with DGW, therefore you would have to upgrade it using pip in order to run it.
```
pip install --upgrade pandas
```

Packages that are not in EPD also need to be installed using pip
```
pip install pysam fastcluster
```
Don't forget MLPY:
```
pip install -e git://github.com/sauliusl/mlpy.git#egg=mlpy
```




