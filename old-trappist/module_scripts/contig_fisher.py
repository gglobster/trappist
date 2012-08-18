__author__ = 'GG'

## PSEUDOCODE FOR CONTIG FISHER SCRIPT

# SETUP

# Import methods
from analysis.dataset_manipulation import seq_subset_load, genome_set_load

# Set up environment and config parameters
## create config objects


# Prepare query
subset, subset_file = seq_subset_load(infile, subset_mode, subset_args)
print 'something for the logger'

# Prepare genomes / contig sets
genomes_set, g_set_file = genome_set_load(input_file, input_prefs)


# Blast the query against the databases
