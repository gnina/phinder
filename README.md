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
To run docking on 


