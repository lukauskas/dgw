Dynamic Genome Warping (DGW)
===============

Installation
============================
`numpy`, `scipy`, `pandas==0.10.1`, `mlpy`, `pysam`, `fastcluster` and `matplotlib` packages.

Enthought Distribution
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

