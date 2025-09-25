#!/usr/bin/env python3
from setuptools import setup, find_packages, Extension
from Cython.Build import cythonize
from Cython.Distutils import build_ext

# Cython extension
extensions = [
    Extension("orfipy.orfipy_core", ["orfipy/orfipy_core.pyx"])
]

setup(
    ext_modules=cythonize(extensions),
    cmdclass={'build_ext': build_ext},
)
