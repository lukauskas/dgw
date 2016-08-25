from setuptools import setup, Extension
import os


try:
    from Cython.Build import cythonize
except ImportError:
    cython_supported = False
else:
    cython_supported = True

import numpy
np_lib = os.path.dirname(numpy.__file__)
np_inc = [os.path.join(np_lib, 'core/include')]

cmdclass = {}
# Compile packages from mlpy distribution
if cython_supported:
    ext_modules = [Extension("dgw._mlpy.dtw",
                           ["mlpy_src/dtw/cdtw.c",
                            "mlpy_src/dtw/dtw.pyx"],
                           include_dirs=np_inc)]

    ext_modules = cythonize(ext_modules)
else:
    ext_modules = [Extension("dgw._mlpy.dtw",
                            ["mlpy_src/dtw/cdtw.c",
                             "mlpy_src/dtw/dtw.c"],
                            include_dirs=np_inc)]

setup(
    name='dgw',
    version='0.1.1',
    packages=['dgw', 'dgw.bin', 'dgw.cli', 'dgw.cluster', 'dgw.data', 'dgw.data.parsers', 'dgw.data.visualisation', 'dgw._mlpy',
              'dgw.dtw', 'dgw.tests.data.parsers', 'dgw.tests.data', 'dgw.tests.dtw', 'dgw.tests', 'dgw'],
    install_requires=['argparse',
                      'numpy>=1.6.1', 'scipy>=0.9.0', 'pandas>=0.10.1', 'pysam>=0.7.4',
                      'fastcluster>=1.1.7'],
    extras_require ={
        'visualisation': ['matplotlib>= 1.1.0',
                          'palettable>=2.1.1',
                          'seaborn>=0.7.1']
    },
    entry_points={
        'console_scripts': [
            'dgw-extract-gene-regions = dgw.bin.extract_gene_regions:main',
            'dgw-overlaps2poi = dgw.bin.overlaps2poi:main',
            'dgw-prototypes2dot = dgw.bin.prototypes2dot:main [visualisation]',
            'dgw-worker = dgw.bin.worker:main'
        ],
        'gui_scripts': [
            'dgw-explorer = dgw.bin.explorer:main [visualisation]',
        ]
    },
    ext_modules=ext_modules,

    url='http://sauliusl.github.com/dgw/',
    license='GPLv3',
    author='Saulius Lukauskas',
    author_email='luksaulius@gmail.com',
    description='Dynamic Genome Warping',
    cmdclass=cmdclass,
)

