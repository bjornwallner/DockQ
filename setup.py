from setuptools import setup, Extension

setup(
    ext_modules=[
        Extension(
            "DockQ.operations", sources=["src/DockQ/operations.pyx"],
        ),
    ],
)

