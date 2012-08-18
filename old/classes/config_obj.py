__author__ = 'GG'

class FisherConfig(object):
    """Store configuration parameters for ContigFisher run."""
    def main(self):
        print "This is a FisherConfig object."

    def __init__(self, data_dir, query_name, query_infile, init_mode):
        self.data_dir = data_dir
        self.query_name = query_name
        self.query_infile = query_infile
        self.init_mode = init_mode


class SorterConfig(object):
    """Store configuration parameters for ConservSorter run."""
    def main(self):
        print "This is a SorterConfig object."

    def __init__(self,data_dir, aln_dir, matchup_mode, chop_mode, aln_mode,
                 clump_mode, simil_set, sort_prefs):
        self.data_dir = data_dir
        self.aln_dir = aln_dir
        self.matchup_mode = matchup_mode
        self.chop_mode = chop_mode
        self.aln_mode = aln_mode
        self.clump_mode = clump_mode
        self.simil_set = simil_set
        self.sort_prefs = sort_prefs
