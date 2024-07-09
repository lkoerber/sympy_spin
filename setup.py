from setuptools import setup, find_packages

with open("README.md", "r") as f:
    long_description = f.read()

setup(
    name="sympy_spin",
    version="0.0.1",
    description="Symbolic computation package for quantum magnetism built on-top of SymPy",
    #package_dir={"": "sympy_spin"},
    packages=find_packages(),
    py_modules=["sympy_spin"],
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/lkoerber/sympy_spin",
    author="lkoerber",
    license="MIT",
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3 :: Only",
        "Topic :: Scientific/Engineering :: Physics",
        "Topic :: Scientific/Engineering :: Mathematics",
    ],
    install_requires=["sympy"],
    python_requires=">=3.10",
)
