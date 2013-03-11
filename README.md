Dynamic Genome Warping (DGW)
===============

# Installation

## Prerequisites
DGW depends on `numpy`, `scipy`, `pandas==0.10.1`, `mlpy`, `pysam`, `fastcluster`, `matplotlib` and modified `mlpy` packages (see below).

### MLPY
DGW depends on modified version of MLPY, available from https://github.com/sauliusl/mlpy
You can install it by
```
pip install cython
pip install -e git://github.com/sauliusl/mlpy.git#egg=mlpy
```
Cython is required for compilation of mlpy from source.
Note that mlpy depends on GSL (with header-files). 
Please look for instructions on how to install it on your platform if the installation fails.

On ubuntu this can be installed by 
```
apt-get install libgsl0-dev
```

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
Note that some packages, e.g. `scipy` depend on dev versions of some of the GNU libraries, 
you will likely to have to install them before proceeding, for instance, on ubuntu you will need to:
```
apt-get install libblas-dev liblapack-dev gfortran
```
Similarly, you can install `numpy` and `scipy` using apt-get (or other package manager) if you like.

Finally install modified MLPY:
```
pip install -e git://github.com/sauliusl/mlpy.git#egg=mlpy
```

### Mac OS X 
Probably the most reliable way of installing `numpy`, `scipy`, `pandas`, `matplotlib` and `pysam` is via MacPorts (https://www.macports.org/)
This can be achieved by

```
port install python27 py27-numpy py27-scipy py27-pandas py27-matplotlib py27-pysam 
```

Then the remaining dependancies, `fastcluster` and `mlpy` can be installed from pip
```
pip install fastcluster
pip install mlpy
```

Alternatively, Enthought Python Distribution can be used to install these packages (see instructions for Windows).
However `matplotlib` may not work correctly as EPD does not install python as a framework

### Windows
It is probably easiest to install `numpy`, `scipy`, `pandas` and `matplotlib` on Windows via Enthought python distribution. This distribution is free for academic use. See 
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

## Installing Python 2.7 to a non-root environment
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

# Installation of DGW

DGW has to be installed directly from the source. You can obtain source by cloning the DGW repository.

```
git clone
```

In order to be able to use DGW from any location in your system
set up your PATH and PYTHONPATH variables to point to the location of the new repository, on unix this is done:
```
export PYTHONPATH=$PYTHONPATH:/directory/where/dgw/is/checked/out
export PATH=$PATH:/directory/where/dgw/is/checked/out/bin
```
Note the "bin" in the end of PATH directory -- this is where the main executables of package lie in. 
You may want to add these two lines to your `~/.bashrc` so they are executed every time you open your shell.

# Usage

DGW is split into two parts - computationally demanding part, `dgw-worker` and an exploratory part - `dgw-explorer`.

## `dgw-worker`

The worker part of the module is responsible for the actual hard work done in clustering the data. It preprocesses the data, computes intermediate representations, calculates DTW distances between the data, performs hierarchical clustering and calculates prototypes of the clusters.

### Sample usage

Typically, `dgw-worker` would be run as follows:
```
dgw-worker.py -r tss_regions.bed  -d wgEncodeBroadHistoneK562H3k4me3StdAlnRep1.bam wgEncodeBroadHistoneK562Pol2bStdAlnRep1.bam --prefix dgw_example
```
In this case we are providing a bed file of regions of interest we want to cluster (`-r tss_regions.bed`), two datasets to work on (`-d wgEncodeBroadHistoneK562H3k4me3StdAlnRep1.bam wgEncodeBroadHistoneK562Pol2bStdAlnRep1.bam`) and setting the prefix of files that will be output to `dgw_example`.

The DGW-worker will take all alignments from both datasets at regions in the `tss_regions.bed`. These alignments will then be extended and put into bins of 50 base pairs wide (use `-res` parameter to change this width).
Then the unexpressed regions that have no bin with more than 10 reads in them (`-min-pileup` constraint to change) will be ignored. Note that these ignored regions are then saved to `{prefix}_filtered_regions.bed` file.
The remaining data will be normalised by adding two artificial reads for each bin and then taking the log of the number of reads in the bins.
The remaining regions will then be clustered hierarchically using DTW distance with default parameters.

### Output
The worker will output 8 files to the working directory where `{prefix}` is the prefix specified by `--prefix` argument.

* `{prefix}_config.dgw` -- The main file storing the configuration of DGW that was used to produce the other files.
* `{prefix}_dataset.pd` -- Processed dataset after the normalisation. This can then be passed in a subsequent DGW session as `--processed-dataset` parameter.
* `{prefix}_filtered_regions.bed` -- Regions that were filtered out of the original regions set due to preprocessing constraints.
* `{prefix}_linkage.npy` -- Precomputed linkage matrix that is used in hierarchical clustering 
* `{prefix}_missing_regions.bed` -- regions that were in the BED file provided as an input, but were not in one of the BAM files.
* `{prefix}_prototypes.pickle` -- computed prototypes of the clusters
* `{prefix}_regions.pd` -- regions that were processed, saved in DGW-readable format
* `{prefix}_warping_paths.pickle` -- computed warping paths of the original data projected onto prototypes

### Points of interest
In some cases one would want to track some points of interest and their locations after warping. `dgw-worker` can also be run with a `-poi` parameter specified, for instance:

```
dgw-worker.py -r tss_regions.bed -poi first_splicing_site_locations.bed  -d wgEncodeBroadHistoneK562H3k4me3StdAlnRep1.bam wgEncodeBroadHistoneK562Pol2bStdAlnRep1.bam --prefix dgw_example
```

The regions in `first_splicing_site_locations.bed` must have the same names as the regions in `tss_regions.bed` otherwise DGW won't be able to match them. Also have a look at `--ignore-poi-non-overlaps` id some of the regions in the input file may not contain some of the regions listed as points of interest.
Similarly, `--ignore-no-poi-regions` will make DGW ignore those regions in input file that do not contain any of the points of interest provided.

## `dgw-explorer`




