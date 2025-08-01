# Evaluation of Microbiome HiFi Amplicon Sequencing Simulator (MHASS)

Evaluation scripts and data for all analyses in our evaluation of MHASS (https://github.com/rhowardstone/MHASS)

## Installation

```bash
git clone https://github.com/rhowardstone/MHASS_evaluation.git
cd MHASS_evaluation
bash scripts/setup.sh
```
Download this data file and unzip it to the same directory:
https://drive.google.com/file/d/1F0I_pwu1e4gKK796CsgUZGGq2kQi3Bs0/view?usp=sharing

Then download the ATCC dataset: https://github.com/PacificBiosciences/DevNet/wiki/16S-Data-Set-Sequel-II-System-2.0-Release

And the Phylotag dataset: https://genome.jgi.doe.gov/portal/PhyloTag/PhyloTag.home.html

After installation, update the references in the runner scripts to reflect the location of your downloads, and run:

```bash
bash scripts/1A_runner.sh
```
