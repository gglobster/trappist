from analysis.blasting import blast_record_set

__author__ = 'GG'

class GenomeSet(object):
    """Store data relative to a whole genome or set of contigs."""
    def main(self):
        print "This is a GenomeSet object."

    def __init__(self, in_file, name):
        self.in_file = in_file
        self.name = name
        self.match_sets = []

    def run_blast(self, blast_run):
        query_segs = blast_run.query.query_segs
        blast_prefs = blast_run.parameters
        matches_dict = blast_record_set(self.name, query_segs, blast_prefs)
        self.match_sets.append({'run_id': blast_run.run_id,
                             'matches': matches_dict})

class FisherQuery(object):
    """Store data relative to a query used to blast a genome set."""
    def main(self):
        print "This is a FisherQuery object."

    def __init__(self, query_id, query_segs, query_file):
        self.query_id = query_id
        self.query_segs = query_segs
        self.query_file = query_file

class BlastRun(object):
    """Store data relative to a blast run."""
    def main(self):
        print "This is a BlastRun object."

    def __init__(self, run_id, query, parameters):
        self.run_id = run_id
        self.query = query
        self.parameters = parameters

class GenomicStrip(object):
    """Store data from a single DNA entity."""
    def main(self):
        print "This is a GenomicStrip object."

    def __init__(self, index, name, seqfile, offset, ordinal):
        self.index = index
        self.name = name
        self.seqfile = seqfile
        self.offset = offset
        self.ordinal = ordinal
        self.seg_sets = []

    def load_seqfile_contents(self):
        # get sequence file contents
        from methods.analytics.sequence_file_processes import ensure_fasta,\
            load_fasta
        self.fasta_f = ensure_fasta(self.seqfile)
        self.seqrecord = load_fasta(self.fasta_f)
        return self.seqrecord

    def add_segments_set(self, segment_set):
        # store segment data from coordinates file
        self.segments   = segment_set
        return self