from setuptools import setup, Extension
from Cython.Build import cythonize
import numpy

extensions = [
    Extension(
        "DockQ.operations",
        ["src/DockQ/operations.pyx"],
        include_dirs=[numpy.get_include()],
    ),
]

setup(
    name="dockq",
    ext_modules=cythonize(extensions),
    package_data={
        "src/DockQ": ["operations.pyx"],
    },
)
