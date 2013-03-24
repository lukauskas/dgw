DGW: Dynamic Genome Warping
===============================

Dynamic Genome Warping (DGW) is  an open source clustering and alignment tool for epigenomic marks.
DGW Utilises Dynamic Time Warping distance to adaptively rescale the matching genomic marks to capture similarities
based on their shapes.

DGW is written in Python and needs Python 2.7 to run.
You can check the version of python you are running by doing `python -V` in your terminal.
If you do not have Python installed, or your python is not 2.7, please install it.
See section below how to do this if you do not have root access in the system.

Dependencies
-------------------------------
DGW depends on a fair share of popular open source packages:

- `numpy` - used for standard numerical tasks,
- `scipy` - used for for its cluster.hierarchy module
- `fastcluster` - used for for an efficient implementation of some of the functions of that very same module.
- `pandas` - used for for efficient data storage and processing containers
- `matplotlib` - used for for visualisation
- `pysam` - used for for SAM file processing

The package also uses a modified version of `mlpy`s DTW package (distributed together with the package).

Installation
===============================
Stable installers of DGW are distributed over PyPi: https://pypi.python.org/pypi/dgw


In most cases `easy_install` should be able to install DGW and its dependancies directly by typing::

    easy_install dgw

If this fails, continue reading the following instructions on how to install DGW from source or in an environment without root access.

Platform-specific instructions
-------------------------------
If the easy installation fails, due to dependancies that were not satisfied try installing them manually using
`easy_install` or `pip`, e.g.::

    pip install numpy
    pip install scipy
    pip install pandas
    pip install pysam
    pip install fastcluster
    pip install matplotlib

and then restart the installation.

Ubuntu
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Note that some packages, e.g. `scipy` depend on dev versions of some of the GNU libraries.
You will likely to have to install them before proceeding, for instance, on ubuntu you will need to perform::

    apt-get install libblas-dev liblapack-dev gfortran

Mac OS X
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
If your system is running Mac OS X, MacPorts can be used to install the majority of these dependencies::

    port install python27 py27-numpy py27-scipy py27-pandas py27-matplotlib py27-pysam

The remaining packages can then be installed with `pip` ::

    pip install fastcluster


.. _MacPorts = https://www.macports.org/

Windows
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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

Installing Python 2.7 to a non-root environment
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The above steps assume you have Python 2.7 installed in your system.
If you do not have Python2.7, you need to install it.
The following steps show how to install python to an environment on linux you do not have root access to.

1. Download and extract Python sources::

    wget http://www.python.org/ftp/python/2.7.3/Python-2.7.3.tgz
    tar xvf Python-2.7.3

2. Install python to a local directory::

   cd Python-2.7.3
   ./configure
   make altinstall prefix=~/python_dev/python/ exec-prefix=~/python_dev/python

where `~/python_dev/python` is the desired location to install python to (change as appropriate).
Note the tilde (`~`) indicating this is under `$HOME` directory -- directory my user has access to.
At this point you should have a `python2.7` executable at `~/python_dev/python/bin/`.

3. Set up PATH variables::

    export PATH=~/python_dev/python/bin:$PATH
    export PYTHONPATH=~/python_dev/python/lib/python2.7/site-packages/

4. Download and install setuptools. Make sure that your `PYTHONPATH` variable is set correctly before doing this::

    wget https://pypi.python.org/packages/2.7/s/setuptools/setuptools-0.6c11-py2.7.egg#md5=fe1f997bc722265116870bc7919059ea
    sh setuptools-0.6c11-py2.7.egg

5. You now should be able to install pip by::
    easy_install-2.7 pip

6. Once pip is installed, you can install DGW as usual.
Make sure you use the newly installed `pip-2.7`, which will be in your local directory and not the one that comes with system.


Installation from source
~~~~~~~~~~~~~~~~~~~~~~~~~~
If you want to get the latest version of DGW, obtain the latest source by cloning the repository::

    git clone git://github.com/sauliusl/dgw.git

Navigate to the newly created `dgw` directory and run the installation script::

    python setup.py install


Usage
=======================

DGW is split into two parts - computationally demanding part, `dgw-worker` and an exploratory part - `dgw-explorer`.

`dgw-worker`
-----------------------

The worker part of the module is responsible for the actual hard work done in clustering the data.
It preprocesses the data, computes intermediate representations, calculates DTW distances between the data,
 performs hierarchical clustering and calculates prototypes of the clusters.

Sample usage
~~~~~~~~~~~~~~~~

Typically, `dgw-worker` would be run as follows:
```
dgw-worker.py -r regions.bed  -d dataset1.bam dataset2.bam --prefix dgw_example
```

In this case we are providing a bed file of regions of interest we want to cluster (`-r regions.bed`),
two datasets to work on (`-d dataset1.bam dataset2.bam`) and setting the prefix of files that will be output to `dgw_example`.

The DGW-worker will take all alignments from both datasets at regions in the `regions.bed`.
These alignments will then be extended and put into bins of 50 base pairs wide (use `-res` parameter to change this).
Then the unexpressed regions that have no bin with more than 10 reads in it (`-min-pileup` constraint to change) will be ignored.
Note that these ignored regions are then saved to `{prefix}_filtered_regions.bed` file.
The remaining data will be normalised by adding two artificial reads for each bin and then taking the log of the number of reads in the bins.
The remaining regions will then be clustered hierarchically using DTW distance with default parameters.

Output
~~~~~~~~~~~~~~~~~~~~~~
The worker will output 8 files to the working directory where `{prefix}` is the prefix specified by `--prefix` argument.

* `{prefix}_config.dgw` -- The main file storing the configuration of DGW that was used to produce the other files.
* `{prefix}_dataset.pd` -- Processed dataset after the normalisation. This can then be passed in a subsequent DGW session as `--processed-dataset` parameter.
* `{prefix}_filtered_regions.bed` -- Regions that were filtered out of the original regions set due to preprocessing constraints.
* `{prefix}_linkage.npy` -- Precomputed linkage matrix that is used in hierarchical clustering
* `{prefix}_missing_regions.bed` -- regions that were in the BED file provided as an input, but were not in one of the BAM files.
* `{prefix}_prototypes.pickle` -- computed prototypes of the clusters
* `{prefix}_regions.pd` -- regions that were processed, saved in DGW-readable format
* `{prefix}_warping_paths.pickle` -- computed warping paths of the original data projected onto prototypes

Points of interest
~~~~~~~~~~~~~~~~~~~~~
In some cases one would want to track some points of interest and their locations after warping,
for instance, we might want to see where transcription start sites are mapped to after the warping.
To do this, `dgw-worker` need to be run with a `-poi` parameter specified, for instance::

    dgw-worker.py -r regions.bed -poi poi.bed  -d dataset1.bam dataset2.bam --prefix dgw_example

The regions in `poi.bed` must have the same names as the regions in `tss_regions.bed` otherwise DGW won't be able to match them.
Also have a look at `--ignore-poi-non-overlaps` id some of the regions in the input file may not contain some of the regions listed as points of interest.
Similarly, `--ignore-no-poi-regions` will make DGW ignore those regions in input file that do not contain any of the points of interest provided.

Runtime
-----------------
Please note that DGW Worker is a very computationally-demanding piece of software.
It is designed to be used on a performant computer with as much CPU cores as possible.

A good way to estimate how long will the computation take on your machine is to use `--random-sample` parameter, e.g. pass `--random-sample 1000`.
This parameter will take only a random sample of N regions, where N is the provided number (in this case 1000).
The DGW worker will work on this random sample and report you both the time it took to compute the pairwise distances
on the random sample, and the estimated time to compute them on the full sample.

Prototype estimation and DTW projections onto prototypes will take around an extra 50% of time taken for pairwise distance calculations.

`dgw-explorer.py`
----------------------

A second major part of DGW is the DGW explorer.
This software is much less computationally demanding than DGW Worker and is designed to allow you to explore the results.

In order to use it start it by passing a `{prefix}_config.dgw` file computed by
DGW worker:
``dgw-explorer.py dgw_config.dgw``

The remaining files output by DGW explorer must be in the same directory as the `dgw_config.dgw` file, otherwise the explorer will not be able to locate them.

Upon successful start, a window showing the dendrogram and heatmap will pop up. Left click on the dendrogram to cut it at the desired place, wait for the plot to refresh and click preview to bring up a cluster explorer.

The cluster explorer allows you to cycle through clusters generated by the dendrogram cut and save both the data of the clusters and the generated heatmaps.

Note that you can also provide `-poi` parameter to `dgw-explorer.py`.
This will override the points of interest specified by worker.
DGW Explorer allows you to specify up to two sets of points of interest (just add the -poi parameter twice).

The resulting points of interest will be plotted on top of the heatmap in cluster viewer.

Using the GUI: Main Window
~~~~~~~~~~~~~~~~~~~~~~
When `dgw-explorer` will first start, you would see a dendrogram generated from hierarchical clustering and the
heatmap of data besides it.

*Right click* anywhere on the dendrogram to cut it. You should see the colours of the dendrogram nodes change
corresponding to different cluster assignments. Please be patient if you are working with large dataset, as it
might take a while to redraw the dendrogram.

Once the dendrogram is drawn, clicking on the Preview button, will bring the cluster explorer window up.

Using the GUI: Cluster Explorer
~~~~~~~~~~~~~~~~~~~~~~~




Utility modules
---------------
DGW comes with a few utility modules to help in experiments.

`dgw-extract-gene-regions`
~~~~~~~~~~~~~~~~~~~~~~~~~~
One of these modules is `dgw-extract-gene-regions`. As the name suggests it allows extraction of regions related to
genes. Currently it supports only knownGene files downloaded from ENCODE_.

To obtain these files, navigate to `table browser`_, make sure group `Genes and Gene Prediction Tracks` is selected,
track is set to `UCSC Genes` and the table set to `knownGenes`.

Click `Get output` to get the data and save it to file.

This file can then be processed using `dgw-extract-gene-regions`.


This utility takes two filenames, one for input file, other for the output file and one of the three options
as input parameters:

   - `--gene` - return regions spanning the length of whole gene
   - `--exon N` - return Nth exon (numbering is zero-based, so 0 is first exon).
   - `--splicing-site N` -- return Nth splicing sites (again 0 based).
   - `--tss` -- return transcription start sites of the genes only.

The last two options, `--splicing-site N` and `--tss` take an optional parameter `--window WINDOW_SIZE` that allows the user
to get the window of `WINDOW_SIZE` base pairs around the data.

For instance, if we wanted to get regions with 2000 bp around the transcription sites of all known genes::

    dgw-extract-gene-regions --tss --window 2000 knownGenes regions_around_tss.bed

The resulting regions will be saved to `regions_around_tss.bed`.

If one wants just the locations of TSS, i.e. for visualising in dgw-explorer, specify a window size of zero: `--window 0`.


.. _ENCODE: http://encodeproject.org/ENCODE/
.. _table browser: http://encodeproject.org/cgi-bin/hgTables?hgsid=330609261&clade=mammal&org=Human&db=hg19&hgta_group=genes&hgta_track=wgEncodeRegTxn&hgta_table=0&hgta_regionType=genome&position=chrX%3A151073054-151383976&hgta_outputType=wigData&hgta_outFileName=

`dgw-overlaps2poi`
~~~~~~~~~~~~~~~~~~~~
DGW also comes with an utility to ease the generation of POI files.
This utility, `dgw-overlaps2poi` takes two bed files for input:

- `main_regions_of_interest.bed`, the file that contains
the regions that will be processed by DGW-Worker (i.e. peak caller results)
- `poi_regions.bed`, bed file containing a list of points of interest, e.g. locations of transcription start sites
generated by dgw-extract-gene-regions.

The utility will then process all the regions in `main_regions_of_interest.bed` find all the regions in `poi_regions.bed`
that overlap *completely* (i.e. *all* points in the points of interest region are contained within the main region)
and spit out a DGW-readable POI file to standard output (which then can be redirected to file).

Example usage::

    dgw-overlaps2poi macs_results.bed tss_regions.bed > tss.poi

`dgw-prototypes2dot`
~~~~~~~~~~~~~~~~~~~~~
This tool is a helper tool around `dgw-explorer` that helps to visualise and debug the prototype generation out of
the DGW result config. After running it will generate a set of images, corresponding to original data and created prototypes
and a dot file, that can then be converted to other formats using GraphViz_.

Please consult::

   dgw-prototypes2dot --help

for more information.

.. WARNING::
   This function will generate png images for every single node in the dendrogram, including the leaves (actual data).
   This means that might take a fair amount of time to run and generate gigabytes of data, therefore it is not recommended
   to run it for anything but small datasets in order to understand how prototype generation works.

.. _GraphViz: http://www.graphviz.org/

Quickstart
=======================
This section will walk you though some example usage of DGW in full so you can have running start with the software
It is going to use MACS peak caller and the dataset