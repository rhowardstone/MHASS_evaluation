#!/bin/bash
# setup.sh
# Downloads and installs MHASS, additional dependencies.


unzip MHASS_datasets.zip

pip install scipy edlib tqdm matplotlib

git clone https://github.com/rhowardstone/MHASS.git
cd MHASS
bash install_dependencies.sh
conda activate mhass
cd ../




