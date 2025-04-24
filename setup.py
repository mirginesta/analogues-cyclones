import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="analogues-cyclones",
    version="0.1.0",
    author="Mireia Ginesta",
    author_email="Mireia.Ginesta@smithschool.ox.ac.uk",
    description="Python package for computing analogues of cyclones",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/mireiaginesta/analogues-cyclones",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
    install_requires=[
        'numpy',
    ],
)
