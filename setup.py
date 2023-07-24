import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="CellNeighborEX", 
    version="0.0.3", 
    author="hkim240", 
    author_email="hbkim20005@gmail.com", 
    description="library for neighbor-dependent gene expression analysis", 
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/hkim240/CellContact", 
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.7',
)