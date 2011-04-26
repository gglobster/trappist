__author__ = 'GG'

class GenomicStrip(object):
    """Store data from single DNA entity (plasmid, genome region etc.)."""
    def main(self):
        print "This is a GenomicStrip object."

    def __init__(self, index, name, seqfile, offset, ordinal):
        self.index = index
        self.name = name
        self.seqfile = seqfile
        self.offset = offset
        self.ordinal = ordinal

    def load_seqfile_contents(self):
        """Get sequence file contents."""
        from methods.analytics.sequence_file_processes import ensure_fasta,load_fasta
        self.fasta_f = ensure_fasta(self.seqfile)
        self.seqrecord = load_fasta(self.fasta_f)
        self.sequence = self.seqrecord.seq
        self.length = len(self.sequence)    # consider throwing this out
        return self.seqrecord

    def add_segments_set(self,segment_set):
        """Store segment data from coordinates file."""
        self.segments   = segment_set
        return self

class GenomeSet(object):
    """Store data relative to a whole genome or set of contigs."""
    def main(self):
        print "This is a GenomeSet object."

    def __init__(self,in_file,name):
        self.in_file = in_file
        self.name = name
        self.contigs = []
        self.sorted = []

        