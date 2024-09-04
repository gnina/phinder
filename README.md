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
Currently, Phinder is meant to be run within a project directory of which it's structure is very important for easily working with various files generated from the program. The best way to run phinder is to do the following:

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
  
Phinder clusters the docked fragments to generate six pharmacophore feature types:
  * Hydrogen bond donor
  * Hydrogen bond acceptor
  * Aromatic
  * Hydrophobic
  * Positive Ion
  * Negative Ion

## Phinding

To phind pharmacophores use the Phinder.py file on your receptor of choice. Phinder will auto-generate the binding pocket from an input ligand. Alternatively, you can provide Cartesian coordinates to pocket. Phinder will dock the fragments and generate a file of all pharmacophore features called GeneratedPharma.json. To run phinder use the following command:

```bash
## From Phinder Home Directory using ligand to auto-generate a pocket. The output directory defaults to the current directory. 
python Phinder.py -r receptor -l ligand -w output_directory
```
```bash
## From Phinder Home Directory using Cartesion coordinates to generate a pocket. 
## You must provide the center X,Y, and Z coordinates as well as the size in each coordinate
python Phinder.py -r receptor -p="coords" --coords=[X-coordinate,Y-coordinate,Z-coordinate,X-size,Y-size,Z-size] -w output_directory
```

## Machine Learning
Phinder includes the option to utilize an MLP classifier trained on the PDB-Bind refined set (http://www.pdbbind.org.cn/index.php).  If this optin is selected, Phinder will automatically output a json of the top 10 features ranked by this model called Top10.json.  
```bash
python Phinder.py -r receptor -l ligand -w output_directory -m True 
```

## Exhaustiveness 
Gnina performs docking through Monte Carlo sampling. The maximum number of retained conformations per fragment defaults to 100, but can be set by the user. 
```bash
#n=25 means a maximum of 25 poses per fragment will be retained
python Phinder.py -r receptor -l ligand -w output_directory -n 25 
```
Congrats! You have completed your first Phind! Now you can review your Pharmacophores in PyMOL. To do so make sure you have the load_query module for PyMOL it can be found [here](https://sourceforge.net/projects/pharmer/files/)
