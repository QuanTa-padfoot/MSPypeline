conda env create -f .\environment.yml

conda activate mspypeline_dev

python .\setup.py install

Write-Output "`nInstallation completed. Please remember to define the R_HOME variable before starting the GUI.`n"
