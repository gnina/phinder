# phinder
Under active development - **do not use**

Semi-automatic identification of pharmacophore features from unbound receptor structure using gnina fragment docking.  

## Installation
Currently installation should be done by cloning this GitHub repository as well as following the installation instructions for GNINA 1.0.

To install simply do the following command:

```bash
git clone https://github.com/gnina/phinder.git
```
Once you have cloned the repository you can work with the phinder scripts to generate pharmacophore!

## Usage
Phinder is separted into two main tasks docking and phinding. Currently, Phinder is meant to be run within a project directory of which it's  structure is very important for easily working with files various files generated from the program. The best way to run phinder is to do the following:

1. Create a project and copy Phinder into said project
```bash
mkdir projectname
cp Phinder projectname/Phinder
```
2. Go into said project
```bash
cd projectname
```

## Docking
Docking is done using GNINA 1.0 where we have locked in a set of probes for the fragment docking. Those probes are as follows:
  * Acetate
  * Acetamide
  * Benzene
  * Guanidine
  * Imidazole
  * Isobutane
  * Isopropanol
  * Isopropylamine
  * Water
  

To run docking on with GNINA on all the probes simple use the docker.py module inside the docking package while specifying the receptor (.pdb file type) you would like to dock on, the ligand (.mol2 file type) you want GNINA to autobox with, and the output directory were you would like to send the Docked Probes (Automatically named probeDocked.sdf). **NOTE: As a generate convention in Phinder individual receptor-ligand pairs are kept inside a directory with the name of the receptor, in the PharmacoFindWorkplace. Thus in the example below you can see the relative file paths to said directory.**

```bash
python docker.py -r ../PharmacoFindWorkplace/receptor/receptor.pdb -l ../PharmacoFindWorkplace/receptor/crystal-ligand.mol2 -w ../PharmacoFindWorkplace/receptor/
```

## Phinding

To phind pharmacophores navigate to the Pharmacophore Workplace and run the POpt2.py file inside your receptor directory. POpt2.py will automatically generate individual pharmacophore feature files (featureGeneratedPharma.json) in the OptimalPharmacophore directory within the receptor directory you created for you receptor of interest.

```bash
## From Phinder Home Directory
cd PharmacoFindWorkplace/receptor
python ../POpt2.py
```

You will see the OptimalPharmacophore directory appear in the receptor directory with the generated pharmacophores!

Congrats! You have completed your first Phind! Now you can review your Pharmacophores in PyMOL. To do so make sure you have the load_json module for PyMOL it can be found [here](https://sourceforge.net/projects/pharmer/files/)
