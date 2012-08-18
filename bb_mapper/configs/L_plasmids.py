# Proximity thresholds for clumping
prox_D = 2000   # for ballpark estimation
prox_F = 500    # for fine alignment (set to 0 to skip clumping)

# Size threshold for drawing shading segment
min_size = 100

# Chopping size to limit length of detailed alignments
max_size = 5000
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

# Mauve executable
mauve_exec = 'progressiveMauve'
#mauve_exec = 'mauveAligner'

# Datatypes
mtype = [('A','i4'), ('B','i4'), ('C', 'i4'), ('D','i4')]
segtype = [('A','i4'), ('B','i4'), ('C', 'i4'), ('D','i4'), ('E','i4')]
