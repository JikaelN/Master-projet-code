import pyslim
import os
import msprime
import tskit
import numpy as np
import random
import sys

random.seed(123)

#tree_files = "./R_analysis/500_rep/raw_trees/stepping/"  # First argument: path to the .trees file
#vcf_output_dir = "./R_analysis/500_rep/vcf/stepping/" # Second argument: output directory for VCF files

tree_file = sys.argv[1]
vcf_output_dir = sys.argv[2]

# Ensure output directory exists
os.makedirs(vcf_output_dir, exist_ok=True)
#for tree_file in os.listdir(tree_file):
#    # Vérifie que c'est bien un fichier .trees
    #if tree_file.endswith(".trees"):
        # Construit le chemin complet du fichier
        #tree_file_path = os.path.join(tree_file, tree_file)
    
        #Charge la séquence d'arbres
tree_filename = os.path.basename(tree_file)
replicate = os.path.splitext(tree_filename)[0].split("_")[1]  # Extract replicate number from filename
try:
    ts = tskit.load(tree_file)
    #print(f"Successfully loaded {tree_file}")
except Exception as e:
    sys.exit(1)

# Recapitation with migration between our population
demography = msprime.Demography.from_tree_sequence(ts)
for pop in demography.populations:
    # Put effective population size
    pop.initial_size = 1000
    
# Create ancestral population
demography.add_population(name="ancestral", initial_size=8000)

# Add split: p1 through p8 came from 'ancestral' at time=12000 (SLiM tick)
demography.add_population_split(
time=2000,
derived=[f"p{i}" for i in range(1, 9)],
ancestral="ancestral"
)

migration_rate = 0.001 # 0.2% migration rate
for i in range(1, 9):  # Populations p1 to p8

    if  i < 8:  # Avoid self-migration (no migration within the same population)
        demography.add_migration_rate_change(
            time=2000,
            rate=migration_rate,
            source=f"p{i}",
            dest=f"p{i+1}"
        )
        demography.add_migration_rate_change(
            time=2000,
            rate=migration_rate,
            source=f"p{i+1}",
            dest=f"p{i}"
        )
        
rts = pyslim.recapitate(
    ts, demography=demography,
    recombination_rate=1e-7,
    random_seed = 4
)
# Simplify to keep only 20 individuals per populations
# Get all individuals that are alive at generation 0
alive_inds = pyslim.individuals_alive_at(rts, 0)

# Convert alive individuals to a dictionary by population
pop_dict = {pop_id: [] for pop_id in range(1,9)}  # Store individuals by population

for ind_id in alive_inds:
    ind = ts.individual(ind_id)
    pop_id = ts.node(ind.nodes[0]).population  # Get population ID from first node
    pop_dict[pop_id].append(ind)

# Sample exactly 20 individuals per population
sampled_nodes = []

for pop_id in range(1,9):  # 8 populations
    inds_in_pop = pop_dict[pop_id]

    if len(inds_in_pop) < 20:
        raise ValueError(f"⚠️ Population {pop_id} has only {len(inds_in_pop)} individuals; at least 20 are required.")

    sampled_inds = np.random.choice(inds_in_pop, 20, replace=False)

    # Add both diploid nodes for each sampled individual
    for ind in sampled_inds:
        sampled_nodes.extend(ind.nodes)  # Add both chromosomes

# Simplify tree sequence with the selected sample
simplified_ts = rts.simplify(samples=sampled_nodes)

#print(f"✅ Simplification done. {len(sampled_nodes)//2} individuals ({len(sampled_nodes)} nodes) retained.")

#Add mutation
next_id = pyslim.next_slim_mutation_id(simplified_ts)
mutated_ts = msprime.sim_mutations(
    simplified_ts,
    rate=1e-7,
    model=msprime.SLiMMutationModel(type=0, next_id=next_id),
    keep=True
)
#print(f"The tree sequence now has {mutated_ts.num_mutations} mutations,\n"
#    f"and mean pairwise nucleotide diversity is {mutated_ts.diversity():0.3e}.")

#save as vcf
nts = pyslim.generate_nucleotides(mutated_ts)
nts = pyslim.convert_alleles(nts)
vcf_basename = os.path.basename(tree_file).replace(".trees", ".vcf")
vcf_filename = os.path.join(vcf_output_dir, vcf_basename)
with open(vcf_filename, "w") as vcf_file:
    nts.write_vcf(vcf_file,position_transform = lambda x: np.fmax(1, x))

# Summary
#print("\n✅ Mutations added successfully.")
#print("Number of mutations:", mutated_ts.num_mutations)
#print("Number of segregating sites:", mutated_ts.num_sites)
