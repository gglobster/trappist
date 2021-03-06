# Project data
project_id = 'BCSL'
project_date = '2012'
prot_db_name = 'Bacteria_prot'

from genome_sets.phage_trio import all as genomes
from genome_sets.references import tectiviri as references

# segment context capture
capture_span = 500

# re-annotation flag -- only set to False is the original ref format is gbk
ref_annot_flag = False

# Blast parameters
blast_prefs = {'evalue': 0.01, 'outfmt_pref': 6}
min_nt_match = 300     # min size for a NT blast hit to be relevant
min_nt_score = 100     # min score (NT)
min_nt_idp = 50        # min identity percentage (NT)
min_aa_match = 10      # min size for an AA blast hit to be relevant
min_aa_score = 50      # min score (AA)
min_aa_idp = 50        # min identity percentage (AA)

# Proximity thresholds for clumping
prox_D = 2000   # for ballpark estimation
prox_F = 300    # for fine alignment (set to 0 to skip clumping)

# Size threshold for drawing shading segment
min_size = 100

# Chopping size to limit length of detailed alignments
max_size = 3000
chop_mode = 'exact_size' # 'maxsize_bisect','maxsize_divisor','count_divisor'

# Identity percentage cutoffs and color coding
idpt = {95: '#444444',     # top similarity class (HexColor('#444444'))
        85: '#777777',     # upper middle class (HexColor('#777777'))
        70: '#BBBBBB',     # lower middle class (HexColor('#BBBBBB'))
        50: '#DDDDDD',     # low similarity class (HexColor('#DDDDDD'))
         0: '#FFFFFF'}     # lower than cutoff (HexColor('#FFFFFF'))

# Thresholds for binning contig sizes
ctg_thresholds = [100, 500, 1000]

# Function categories and legend (keys MUST be lowercase)

fct_flags = {'mge': ('transposase', 'transposon', 'intron',
                     'tyrosine recombinase', 'dna-invertase',
                     'reverse transcriptase'),
             'rep': ('something else', 'replication'),
             'syn': ('synthase', 'something else'),
             'tox': ('protective antigen', 'lethal factor',
                     'virulence factor'),
             'ger': ('spore germination protein', 'something else'),
             'tra': ('type iv', 'type ii/iv', 'topoisomerase',
                     'dna translocase ftsk', 'conjugal transfer',
                     'conjugation protein'),
             'ctl': ('transcriptional regulator', 'regulatory protein',
                     'regulator', 'transcriptional repressor'
                     'response regulator aspartate phosphatase'),
             'unk': ('uncharacterized', 'conserved domain protein'),
             'def': ('no match', 'hypothetical protein', 'pXO1-')} # default

fct_colors = {'mge': ('#66CC00', 'MGE'),
              'rep': ('#FF9900', 'Replication'),
              'syn': ('#FF00CC', 'Synthesis'),
              'tox': ('#FF0000', 'Pathogenesis'),
              'ger': ('#993333', 'Germination'),
              'tra': ('#6666FF', 'Transfer'),
              'ctl': ('#FFCC00', 'Control'),
              'unk': ('#666666', 'Uncharacterized'),
              'oth': ('#CCCCCC', 'Other'),
              'def': ('#FFFFFF', 'No match')}

# Concat contigs separator
separator = "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"

