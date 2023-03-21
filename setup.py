from setuptools import setup, find_packages

with open("README.md", "r", encoding = 'utf-8') as fh:
    long_description = fh.read()

with open("mspypeline/version.py", "r", encoding = 'utf-8') as f:
    version = f.readline().split()[-1].strip('"')

setup(
    name="mspypeline",
    version=version,
    description="Package to analyze Mass Spec Data",
    long_description=long_description,
    long_description_content_type="text/markdown",
    packages=find_packages(),
    include_package_data=True,  # include files specified in MANIFEST, i.e. config/
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    install_requires=[
        "numpy>=1.24.0",
        "pandas>=1.5.3",
        "scipy>=1.10.1",
        "rpy2==3.4.5",
        "tzlocal>=4.3",
        #"ruamel_yaml>=0.15.46",
        "matplotlib>=3.7.0",
        "matplotlib-venn>=0.11.9",
        "adjusttext>=0.7.3",
        "scikit-learn>=1.2.2",
    ],
    project_urls={
        "Documentation": "https://mspypeline.readthedocs.io/en/stable/",
        "Source": "https://github.com/siheming/mspypeline",
        "Bug Tracker": "https://github.com/siheming/mspypeline/issues",
    }
)
