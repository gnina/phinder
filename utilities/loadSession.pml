run ~/Downloads/load_query_1_1\ 2/load_query.py
run ~/Downloads/loadDir.py
set all_states, on
load crystal_ligand.mol2
load_query pharma.json
load receptor.pdb
loadDir(suff='.sdf')
load_query Aromaticgeneratedpharma.json
load_query NegativeIongeneratedpharma.json
load_query PositiveIongeneratedpharma.json
load_query Hydrophobicgeneratedpharma.json
load_query HydrogenDonorgeneratedpharma.json
load_query HydrogenAcceptorgeneratedpharma.json
