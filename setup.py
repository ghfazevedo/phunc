from setuptools import setup, find_packages

setup(
    name="phunc",
    version="1.1.0",
    author="Guilherme Azevedo",
    description="A CLI tool to calculate the probability of fixation of differences in a hypothetical nuclear locus that controls phenotype under neutral divergence.",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/ghfazevedo/phunc",
    packages=["phunc"],
    package_dir={"phunc": "src/"},  
    include_package_data=True,
    install_requires=[
        'dendropy',
        'matplotlib',
        'pandas',
        'numpy',
        'scipy'
    ],
    entry_points={
        'console_scripts': [
            'phunc = phunc.phunc:main',
        ]
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: Unix",
    ],
    python_requires='>=3.6',
)

