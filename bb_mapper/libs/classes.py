import sys
from converters import fas2gbk, gbk2fas
from loaders import load_fasta

class Noodle(object):
    """Object that holds info about a piece of DNA."""

    def __init__(self, genome, seq_dir):
        self.name = genome['name']
        self.fas = None
        self.gbk = None
        self.offset = genome['offset']
        self.nudge = genome['nudge']+1
        self.invert = False

        if genome['input'] == 'fas':
            self.fas = seq_dir+genome['file']
            self.gbk = fas2gbk(seq_dir+genome['file'])
        elif genome['input'] == 'gbk':
            self.gbk = seq_dir+genome['file']
            self.fas = gbk2fas(seq_dir+genome['file'])
        else:
            print "ERROR in input format: FastA or Genbank required"
            sys.exit()

        self.len = len(load_fasta(self.fas).seq)
