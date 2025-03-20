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
- Python 3.7, 3.8, or 3.9
- R >= 4.3.0

Installing MSPypeline can be done from an Anaconda Powershell Prompt. First, you should to go to a directory at which you wish to download MSPypeline to (e.g.):
```
cd C:\Temp
```

Running the following lines would download MSPypeline from this repository, set up the virtual environment, as well as install all dependencies:

```
if (Test-Path .\mspypeline) {
     Write-Host "Removing the mspypeline folder from the previous installation."
     Remove-Item .\mspypeline
 }
 
Invoke-WebRequest 'https://github.com/QuanTa-padfoot/MSPypeline/archive/refs/heads/main.zip' -OutFile .\mspypeline.zip

Expand-Archive .\mspypeline.zip –DestinationPath .\
 
Remove-Item .\mspypeline.zip

Rename-Item “MSPypeline-main” –NewName “mspypeline”

cd .\mspypeline

.\setup.ps1
```

*Note: The above commands do not work for some directories (such as C:) due to the operating system, thereby causing errors. In that case, I'd recommend changing the directory and do the setup above again. I've found `C:\Temp` to always work.
After successful installation, the Powershell prompt will print:
```
Please restart your computer.
```

Note: A timeout error sometimes occur due to proxy servers. Please see the instruction from Anaconda [(link)](https://docs.anaconda.com/working-with-conda/reference/security/) to solve this issue.

# Getting started
Before using MSPypeline, you need to define the R file's location in the powershell prompt (e.g.):
```
$env:R_HOME = 'C:\Program Files\R\R-4.3.0'
```

If you have a different version of R, simply replace `4.3.0` with your version.

Now you can start the GUI:
```
conda activate mspypeline_dev
python -m mspypeline --gui
```

# Guide for data analysis
## Accepted data type
MSPypeline supports the analysis of DDA or DIA proteomics data AFTER peptide/protein identification. It accepts DDA data from MaxQuant and DIA data from Spectronaut.

For MaxQuant, all files to analyze need to be in a folder labeled 'txt'.

For Spectronaut, currently accepted file formats are '.xls', '.csv', and '.tsv', though I'd recommend using csv or tsv since it keeps the format of the data better. I'd suggest storing each file in a separate folder (i.e. 1 csv file per folder). The folder would later be used to store data analysis results, hence having several data files might cause confusion as to which file was used for the analysis. The csv file and folder's names do not affect data analysis, so you could name them in a way that is informative and convenient for you.

Data files often use commas as decimal points and/or thousand separators, both of which will be recognized and processed by MSPypeline. Thus, you do not need to manually convert them.

After selecting the directory to the data, select the type of reader depending on the data format: "mqreader" for MaxQuant output, and "SpectroReader" for Spectronaut output.

Most analyses are supported for both data types, except for the MaxQuant report and the peptide analysis.

### Naming convention
Samples should be named as follow: ***Group_Biological Replicate_Technical Replicate***. I'd recommend not having comma, tab, and semi-colon in the names since MSPypeline might wrongly interpret them as separator (and therefore induce an error). Levels should be separated by the "_" sign, and more levels can be added depending on your experiment. All samples must have the same number of levels. 

Avoid having spaces in sample names, as this sometimes causes issue, though the protein and peptide file readers typically remove these spaces when importing files. 

## Renaming samples
Changing the order of the levels in sample names is sometimes desired (e.g. A_C_B instead of A_B_C). To this, click the ***Reorder level*** button and type the name of the first sample in the order you want. After pressing ***OK***, the level change will be applied to all other samples and a file named ***sample_mapping.txt*** will be generated in the config folder to record the old names (i.e. names as in the data file) and new names. Change the "Yaml file" to "default" and press ***Update*** with the correct reader. Double-check if the sample names have changed to your will in the ***sample names*** box.

Please notice that columns in the data file itself will not be changed, only the sample names are changed after mspypeline read the file. 

Alternatively, if you would like to change the name of some samples (e.g. due to typos), you can make a ***sample_mapping.txt*** file in the config folder to instruct MSPypeline how to change the names. The file should have two tab-separated column, the first indicating the name as in the data file and the second indicating the name that you would like to change. For clarity, make two columns containing all sample names in your data (tip: you can take this from the ***config.yml*** file) and then make the change at the target sample(s). 

## MaxQuant report
For quality control of the data, generated via the ***Create report*** button in the lower left corner of the GUI. Only available for MaxQuant data.

## Peptide analysis
Only available for Spectronaut data.

The peptide data should be a .csv or .tsv file and stored in a folder named "peptide" in the data directory. Samples in the peptide data should have identical names as in the protein data. If you have made some changes in the sample names of the protein data file, please make similar changes to the peptide data.

The peptide data file records precursor ions intensity. Peptide intensities are computed by summing intensities of precursors with the same sequence.

For the peptide report, select the protein list(s) from the "Pathway analysis" section. A report will be generated per protein in the selected list and stored in the "peptide_report" folder. Internet connection is required since the protein sequence will be downloaded from UniProt using UniProt ID. A ***ConnectionError*** may occur if there is a proxy issue.

# Bug reporting
This repository is under development, hence errors are likely to occur. Should an error occur, please report it to me at ***congquan.ta@dkfz-heidelberg.de***. In your email, please state:
- What the error is
- A step-by-step procedure leading to the error
- (If applicable) the data file with which you encountered the error
- The error message on the pop-up window and the Anaconda Powershell Prompt. Please also mention if no error message occur, since this is informative!

If you have not updated MSPypeline in a while, I'd recommend re-installing it and see if the error persists. There is a chance that the error is solved in a later version of MSPypeline.


Have fun doing mass spec :D

Quan
