from setuptools import setup, Extension
from Cython.Build import cythonize
import numpy

extensions = [
        Extension(
            "DockQ.operations", sources=["src/DockQ/operations.pyx"],
            include_dirs=[numpy.get_include()],
        ),
    ]

setup(
    ext_modules=cythonize(extensions)
)
