Invoke-WebRequest 'https://github.com/QuanTa-padfoot/MSPypeline/archive/refs/heads/main.zip' -OutFile .\mspypeline.zip

Expand-Archive .\mspypeline.zip –DestinationPath .\
 
Remove-Item .\mspypeline.zip

Rename-Item “MSPypeline-main” –NewName “mspypeline”

cd .\mspypeline

conda env create -f .\environment.yml

conda activate mspypeline_dev

python .\setup.py install

$env:R_HOME = 'C:\Program Files\R\R-4.2.1'

Write-Output "`nPlease restart your computer`n"
