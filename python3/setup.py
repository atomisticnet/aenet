import os
import glob
from setuptools import setup, find_packages, Extension
from Cython.Distutils import build_ext
from numpy import get_include

with open("README.rst") as f:
    long_desc = f.read()
    ind = long_desc.find("\n")
    long_desc = long_desc[ind + 1:]

ext_symmfunc = Extension(
   "aenet.core",
    [os.path.join("aenet", "core.pyx")],
    extra_compile_args=[
        '-I' + os.path.join(os.path.pardir, 'src'),
        '-I' + os.path.join(os.path.curdir, 'aenet'),
        '-fPIC', '-O2'],
    extra_objects = \
       glob.glob(os.path.join(os.path.pardir,"lib","Lbfgsb.*","*_pic.o")) \
     + glob.glob(os.path.join(os.path.pardir,"src","*_pic.o")),
    libraries = ['lapack', 'blas', 'gfortran']
)

setup(
    cmdclass={'build_ext': build_ext},
    name="aenet",
    packages=find_packages(),
    version="0.1.0a1",
    install_requires=["numpy>=1.5"],
    author="Nongnuch Artrith, Alexander Urban",
    author_email="nartrith@atomistic.net, aurban@atomistic.net",
    maintainer="Nongnuch Artrith, Alexander Urban",
    url="http://aenet.atomistic.net",
    license="GNU GPL",
    description="Artificial Neural Network Potentials",
    long_description=long_desc,
    keywords=[],
    classifiers=[
        "Programming Language :: Python :: 2.7",
        "Development Status :: 5 - Production/Stable",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: GNU GPL",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Physics",
        "Topic :: Scientific/Engineering :: Chemistry"
    ],
    include_dirs=[get_include()],
    ext_modules=[ext_symmfunc],
    scripts=glob.glob(os.path.join("scripts", "*.py"))
)
