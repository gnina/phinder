![pharmacophore](https://drive.google.com/uc?export=view&id=1Vk1gMmTzlWvmC2N-DBHXpEqfq_wOkrk9)

# phinder
Under active development - **do not use**

Semi-automatic identification of pharmacophore features from unbound receptor structure using gnina fragment docking.  

## Installation
Currently installation should be done by cloning this GitHub repository as well as following the installation instructions for GNINA 1.0 and rd-kit.

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
cp -r Phinder projectname/Phinder
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
  
Docking is performed if no fragment docked files are detected (fragmentDocked.sdf).

## Phinding

To phind pharmacophores navigate to the Pharmacophore Workplace and run the Phinder.py file inside your receptor directory. Phinder.py will automatically dock and generate a single pharmacophore feature file (GeneratedPharma.json) in the output directory that was specified within the receptor directory you created for you receptor of interest. To run phinder use the following command:

```bash
## From Phinder Home Directory
python Phinder.py -r receptor -l ligand -o output_directory
```

You will see the OptimalPharmacophore directory appear in the receptor directory with the generated pharmacophores!

Congrats! You have completed your first Phind! Now you can review your Pharmacophores in PyMOL. To do so make sure you have the load_query module for PyMOL it can be found [here](https://sourceforge.net/projects/pharmer/files/)
