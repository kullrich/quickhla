# cython: profile=False
# cython: language_level=3

from setuptools import setup
from Cython.Build import cythonize
import numpy as np

setup(
    name="quickhla",
    description="classifies HLA genes",
    version="0.0.1",
    packages=["quickhla"],
    license="MIT"
)
