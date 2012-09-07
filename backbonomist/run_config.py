import numpy

from configs.BCSL import *

# Project root directory
base_root = 'data/'+project_id
g_root_dir = base_root+'/genomes/'
r_root_dir = base_root+'/backbonomist/'

# run-independent directories
fixed_dirs = {
'ref_dbs_dir': base_root+'/ref_dbs/',
'ori_g_dir': g_root_dir+'plasmids/',
'gbk_contigs_dir': g_root_dir+'genbank_ctgs/',
'fas_contigs_dir': g_root_dir+'fasta_ctgs/',
'mfas_contigs_dir': g_root_dir+'mfasta_ctgs/',
'blast_db_dir': g_root_dir+'blast_db/',
'annot_trn_dir': g_root_dir+'annotation/training/',
'ctg_cds_dir': g_root_dir+'annotation/genes/',
'ctg_prot_dir': g_root_dir+'annotation/proteins/',
'ctg_blast_dir': g_root_dir+'annotation/blastp/',
'ctg_stats': g_root_dir+'contig_stats/'
}

# run-dependent directories
run_dirs = {
'pickles': 'pickles/',
'ref_seg_dir': 'references/segments/',
'ref_gbk_dir': 'references/genbank/',
'ref_fas_dir': 'references/fasta/',
'ref_map_dir': 'maps/references/',
'match_pickles': 'matching/',
'blast_out_dir': 'matching/blastn/',
'match_out_dir': 'matching/matches/',
'capture_dir': 'matching/capture/',
'scaffolds_dir': 'scaffolds/',
'run_gbk_ctgs_dir': 'gbk_contigs/',
'mauve_out_dir': 'alignments/mauve_out/',
'aln_seg_dir': 'alignments/aln_segments/',
'maps_dir': 'maps/',
'reports': 'reports/'
}

# Blast results arrays datatypes
blast_dtypes = numpy.dtype([('query', 'S16'),
                           ('dbhit', 'S32'),
                           ('idp', 'float'),
                           ('mlen', 'uint8'),
                           ('mms', 'uint8'),
                           ('gaps', 'uint8'),
                           ('q_start', 'uint16'),
                           ('q_end', 'uint16'),
                           ('r_start', 'uint16'),
                           ('r_end', 'uint16'),
                           ('evalue', 'S5'),
                           ('bitscore', 'float')])

# Mauve executable
mauve_exec = 'progressiveMauve'

# Datatypes
mtype = [('A','i4'), ('B','i4'), ('C', 'i4'), ('D','i4')]
segtype = [('A','i4'), ('B','i4'), ('C', 'i4'), ('D','i4'), ('E','i4')]
