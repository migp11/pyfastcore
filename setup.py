import setuptools
from setuptools import find_packages

install_requires = ['cobra>=0.21.0', 'sympy>=1.0']

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="pyfastcore",
    version="0.0.2",
    author="Miguel Ponce de Leon",
    author_email="miguel.ponce@gbsc.es",
    description="A python-based implementation for the context-specific metabolic model extraction methods from Vlassis et al. 2014",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/migp11/pyfastcore",
    packages=find_packages(exclude=['examples', 'docs', 'tests']),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
    install_requires=install_requires
)
