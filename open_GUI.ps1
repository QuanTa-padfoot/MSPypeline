conda activate mspypeline_dev

$env:R_HOME = 'C:\Program Files\R\R-4.2.1'

python .\setup.py install

python -m mspypeline --gui
