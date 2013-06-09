from setuptools import setup, Extension
import os


try:
    # One needs Cython and numpy to be able to compile the cython code
    from Cython.Build import cythonize
    import numpy
    np_lib = os.path.dirname(numpy.__file__)
    np_inc = [os.path.join(np_lib, 'core/include')]
except ImportError:
    cython_supported = False
else:
    cython_supported = True


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
                             "mlpy_src/dtw/dtw.c"])]

setup(
    name='dgw',
    version='0.1.0',
    packages=['dgw', 'dgw.cli', 'dgw.cluster', 'dgw.data', 'dgw.data.parsers', 'dgw.data.visualisation', 'dgw._mlpy',
              'dgw.dtw', 'dgw.tests.data.parsers', 'dgw.tests.data', 'dgw.tests.dtw', 'dgw.tests', 'dgw'],
    install_requires=['numpy>=1.6.1', 'scipy>=0.10.1', 'pandas>=0.10.1',
              'pysam>=0.7.4', 'fastcluster>=1.1.7', 'matplotlib>=1.2.0'
             ],
    scripts=['bin/dgw-explorer', 'bin/dgw-extract-gene-regions', 'bin/dgw-overlaps2poi',
             'bin/dgw-prototypes2dot', 'bin/dgw-worker'],
    ext_modules=ext_modules,
    url='http://sauliusl.github.com/dgw/',
    license='GPLv3',
    author='Saulius Lukauskas',
    author_email='luksaulius@gmail.com',
    description='Dynamic Genome Warping',
    cmdclass=cmdclass,
)

