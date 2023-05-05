from setuptools import setup
from Cython.Build import cythonize

setup(
    name='Cutils',
    ext_modules=cythonize("cutils.pyx"),
    zip_safe=False,
)