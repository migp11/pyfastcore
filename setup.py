import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="fastcore", # Replace with your own username
    version="0.0.1",
    author="migp11",
    author_email="miguel.ponce@gbsc.es",
    description="A python implementataio fo the context-specific metabolic model extraction ethods com Vlassis et al. 2014",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/migp11/fastcore",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)
