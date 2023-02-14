![full test run](https://github.com/siheming/mspypeline/workflows/full%20test%20run/badge.svg?branch=master)
[![Coverage](https://codecov.io/gh/siheming/mspypeline/branch/master/graph/badge.svg?flag=full-test-run)](https://codecov.io/gh/siheming/mspypeline/branch/master)

![basic test run](https://github.com/siheming/mspypeline/workflows/basic%20test%20run/badge.svg?branch=develop)
[![Coverage](https://codecov.io/gh/siheming/mspypeline/branch/develop/graph/badge.svg?flag=basic-test-run)](https://codecov.io/gh/siheming/mspypeline/branch/develop)

[![Documentation Status](https://readthedocs.org/projects/mspypeline/badge/?version=latest)](https://mspypeline.readthedocs.io/en/latest/?badge=latest)

# README
This pipeline can be used to analyze the results of a MaxQuant analysis.  
To use this package please refer to [the docs](https://mspypeline.readthedocs.io/en/latest/index.html).  

The sources can be found at [PyPI](https://pypi.org/project/mspypeline/) or
[conda](https://anaconda.org/siheming/mspypeline).

# Contribute
If you want to contribute to the project please create a pull request on the develop branch.

# Support
If you encounter a bug or want to request a feature open an issue on github.

# Installation
To use this version of MSPypeline, you would need:
- Anaconda 3
- Python 3.7 or 3.8
- R >= 4.2.1 (for several types of plots)

Installing MSPypeline can be done from an Anaconda Powershell Prompt. First, you should to go to a directory at which you wish to download MSPypeline to (e.g.):
```
cd C:\Temp
```

Running the following lines would download MSPypeline from this repository, set up the virtual environment, as well as install all dependencies:

```
Invoke-WebRequest 'https://github.com/QuanTa-padfoot/MSPypeline/archive/refs/heads/main.zip' -OutFile .\mspypeline.zip

Expand-Archive .\mspypeline.zip –DestinationPath .\
 
Remove-Item .\mspypeline.zip

Rename-Item “MSPypeline-main” –NewName “mspypeline”

cd .\mspypeline

.\setup.ps1
```

After successful installation, the Powershell prompt will print:
```
Please restart your computer.
```

Note: A timeout error sometimes occur due to Internet connection and/or the computer running too many tasks. Re-connecting to the Internet and closing unnecessary files would fix the issue.

# Getting started
Before using MSPypeline, you need to define the R file's location in the powershell prompt (e.g.):
```
$env:R_HOME = 'C:\Program Files\R\R-4.2.1'
```

Now you can start the GUI:
```
conda activate mspypeline_dev
python -m mspypeline –-gui
```
