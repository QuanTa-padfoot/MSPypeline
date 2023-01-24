conda env create -f .\environment.yml

conda activate mspypeline_0.12.0

python .\setup.py install

$env:R_HOME = 'C:\Program Files\R\R-4.2.1'

Write-Output "`nPlease restart your computer`n"