import setuptools
from setuptools import find_packages


with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="pyfastcore",
    version="0.0.6",
    author = "Miguel Ponce-de-Leon",
    author_email = "miguelponcedeleon@gmail.com",
    maintainer = "Miguel Ponce-de-Leon",
    maintainer_email = "miguelponcedeleon@gmail.com",
    description="A python-based implementation for the context-specific metabolic model extraction methods from Vlassis et al. 2014",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/migp11/pyfastcore",
        classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    packages=setuptools.find_packages(where="src"   ),
    package_dir={"": "src"},
    python_requires='>=3.6',
    install_requires=['numpy', 'cobra>=0.25.0']
)
