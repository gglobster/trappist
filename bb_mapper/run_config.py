
from configs.BCSL import *
from configs.L_plasmids import *

# Project root directory
base_root = 'data/'+project_id
g_root_dir = base_root+'/genomes/original/'
r_root_dir = base_root+'/bb_mapper/'

dirs = {
    'seqfiles': g_root_dir,
    'root': r_root_dir,
    'mauve': '/mauve/',
    'aln_segs': '/seg_aln/',
    'maps': '/maps/'
}

