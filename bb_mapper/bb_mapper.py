import sys, time

from libs.writers import ensure_dir
from libs.scripts import align_pairwise, align_multi

from run_config import *
from run_sets import genomes

i_msg = "\n", \
        "##################################################\n", \
        "### BB-Mapper v. 0.2  (prev. PlasmiG)          ###\n", \
        "### Copyright 2012 Geraldine A. Van der Auwera ###\n", \
        "##################################################\n"

def usage():
    print "Usage: python gobble.py bb_map -args (don't call this directly!)"

def info():
    print i_msg

def get_arg(args, key):
    return args[args.index(key)+1]

def main(args=None):

    info()
    
    if args is None:
        args = sys.argv[1:]

    if "-h" in args or "--help" in args:
        usage()
        sys.exit(2)

    if "-id" in args:
        run_id = get_arg(args, "-id")
    else:
        run_id = str(int(time.time())) # use timestamp as unique run identifier

    # determine mode (new alignment or redrawing)
    if "-redraw" in args:
        new_align = False
    else:
        new_align = True

    ensure_dir([dirs['seqfiles'], dirs['root']])
    for item in dirs.values():
        ensure_dir([dirs['root']+run_id+item])

    # one-pair or multiple?
    if len(genomes) < 2:
        raise Exception("ERROR: Need at least two genomes to align!")
    elif len(genomes) == 2:
        align_pairwise(genomes, new_align, dirs, run_id, max_size, chop_mode, mauve_exec)
    else:
        align_multi(genomes, new_align, dirs, run_id, max_size, chop_mode, mauve_exec)

