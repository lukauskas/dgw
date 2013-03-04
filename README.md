Dynamic Genome Warping (DGW)
===============

# Prerequisites
DGW depends on `numpy`, `scipy`, `pandas==0.10.1`, `mlpy`, `pysam`, `fastcluster` and `matplotlib` packages.

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
pip install numpy 
pip install scipy 
pip install pandas
pip install pysam 
pip install fastcluster 
pip install matplotlib
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

# Installing Python 2.7 to a non-root environment
The above steps need python to be installed on the system.
If you do not have Python2.7, you need to install it.
The following steps show how to install python to an environment on linux you do not have root access to.

1. Download and extract Python sources:
```
$ wget http://www.python.org/ftp/python/2.7.3/Python-2.7.3.tgz
$ tar xvf Python-2.7.3
```
2. Install python to local location
```
$ cd Python-2.7.3
$ ./configure
$ make altinstall prefix=~/python_dev/python/ exec-prefix=~/python_dev/python
```
where `~/python_dev/python` is the desired location to install python to (change as appropriate).
Note the tilde (`~`) indicating this is under `$HOME` directory -- directory my user has access to.

At this point you should have a `python2.7` executable at `~/python_dev/python/bin/`.

3. Set up PATH variables
```
$ export PATH=~/python_dev/python/bin:$PATH
$ export PYTHONPATH=~/python_dev/python/lib/python2.7/site-packages/
```

4. Download and install setuptools. Make sure that your `PYTHONPATH` variable is set correctly before doing this.
```
wget https://pypi.python.org/packages/2.7/s/setuptools/setuptools-0.6c11-py2.7.egg#md5=fe1f997bc722265116870bc7919059ea
sh setuptools-0.6c11-py2.7.egg
```

5. You now should be able to install pip by:
```
easy_install-2.7 pip
```

6. Once pip is installed, you can install dependencies for DGW as per linux installations step.
Make sure you use the newly installed `pip-2.7`, which will be in your local directory and not the one that comes with system

# Installation

After the dependencies have been set up, clone this repository:

```
git clone 


