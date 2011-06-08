from analysis.blasting import blast_record_set
from analysis.other_logic import create_id

__author__ = 'GG'

class InfoNote(object):

    """Store optional information about any data object.

    >>> example_note = InfoNote('This is a note.')
    >>> example_note.text
    'This is a note.'

    The note text can be revised with the revise() method.

    >>> example_note.revise('This is a revised version of the note.')
    >>> example_note.text
    'This is a revised version of the note.'

    When a note is revised, the previous version is archived with its
    original timestamp information attached. There is no method provided to
    delete old versions. This is done for provenance tracking purposes.
    
    >>> example_note.versions[0]['text']
    'This is a note.'

    """

    def main(self):
        print "This is an InfoNote object."

    def __init__(self, note_text):
        self.note_id = create_id('info')
        self.text = note_text
        self.versions = []

    def revise(self, new_text):
        self.versions.append({'id': self.note_id, 'text': self.text})
        self.note_id = create_id('info')
        self.text = new_text


class ContigSet(object):

    """Store data relative to a set of genome contigs.

    The basic initialization only stores a link to the input files. This is
    intended to avoid loading the genome data until they are actually needed
    for the analysis.

    >>> contig_set = ContigSet('input_file.fas', 'multifasta')
    >>> contig_set.input_files
    {'multifasta': 'input_file.fas'}

    A BLAST database of the contig set can be created using the make_db()
    method. If the contig sequences are not originally stored in a multiple
    FASTA sequence file, one will be generated automatically using the
    converts2mfasta() method. This method can also be used independently of
    make_db() if necessary for other purposes.

    """

    def main(self):
        print "This is a ContigSet object."

    def __init__(self, input_file, input_type):
        self.set_id = create_id('cs')
        self.input_files = {input_type: input_file}
        self.match_sets = []
        self.notes = []

    def add_name(self, name):
        self.name = name

    def convert2mfasta(self):
        if self.input_files.has_key('multifasta'):
            pass
        elif self.input_files.has_key('multigbk'):
            from Bio import SeqIO
            gbk_file = self.input_files['multigbk']
            seq_file = gbk_file+'.fas'
            SeqIO.convert(gbk_file, 'genbank', seq_file, 'fasta',
                          alphabet=None)
            self.input_files['multifasta'] = seq_file
        else:
            available_formats = self.input_files.keys()
            raise Exception("cannot use available input formats ("+
                             ", ".join(available_formats)+ ")")

    def make_db(self, db_path):
        from analysis.blasting import make_blastDB
        try: self.convert2mfasta()
        except: raise Exception("input file conversion failed")
        else:
            seq_file = self.input_files['multifasta']
        try: make_blastDB(db_path, self.name, seq_file, 'nucl')
        except: raise Exception
        else: self.blast_db = {'db_path': db_path}

    def run_blast(self, blast_run):
        query_segs = blast_run.query.query_segs
        blast_prefs = blast_run.parameters
        try: matches_dict = blast_record_set(self.blast_db, query_segs,
                                             blast_prefs)
        except: raise Exception
        else:
            self.match_sets.append({'run_id': blast_run.run_id,
                                    'matches': matches_dict})

    def add_note(self, note):
        self.notes.append(note)


class BlastQuery(object):
    """Store data relative to a query used to blast a genome set."""
    def main(self):
        print "This is a BlastQuery object."

    def __init__(self, query_id, query_segs, query_file):
        self.query_id = query_id
        self.query_segs = query_segs
        self.query_file = query_file
        self.notes = []

    def add_note(self, note):
        self.notes.append(note)

class BlastRun(object):
    """Store data relative to a blast run."""
    def main(self):
        print "This is a BlastRun object."

    def __init__(self, query, parameters):
        self.run_id = create_id('br')
        self.query = query
        self.parameters = parameters
        self.notes = []
        
    def add_note(self, note):
        self.notes.append(note)

class Segment(object):

    """Store coordinates of a segment of interest.

    >>> segment_1 = Segment(2,5)
    >>> segment_2 = Segment(7,12)
    >>> segment_3 = Segment(13,21)

    """

    def main(self):
        print "This is a Segment object."

    def __init__(self, coord1, coord2):
        self.coord1 = coord1
        self.coord2 = coord2
        self.notes = []

    def add_note(self, note):
        self.notes.append(note)


class SegmentSet(object):

    """Store a list of segments.

    >>> segment_1 = Segment(2,5)
    >>> segment_2 = Segment(7,12)
    >>> segment_3 = Segment(13,21)
    >>> segment_list = [segment_1, segment_2, segment_3]
    >>> segment_set = SegmentSet(segment_list)

    The identifier of the sequence record to which the segments refer can be
    added as source_id for later retrieval in case the segment list is used
    separately from the sequence record.

    >>> from Bio.Seq import Seq, SeqRecord
    >>> example_seq = Seq('AATTTAATATTGCAGATAGGCAGG')
    >>> example_seqrecord = SeqRecord(example_seq, id='example')
    >>> segment_set.add_source_id(example_seqrecord.id)
    >>> segment_set.source_id
    'example'

    """

    def main(self):
        print "This is a SegSet object."

    def __init__(self, segment_list):
        self.segments = segment_list
        self.notes = []

    def add_source_id(self, seqrecord_id):
        self.source_id = seqrecord_id

    def add_note(self, note):
        self.notes.append(note)



class GenomicStrip(object):

    """Store data from a single DNA entity.

    This object is created from a single sequence record object.

    >>> from Bio.Seq import Seq, SeqRecord
    >>> example_seq = Seq('AATTTAATATTGCAGATAGGCAGG')
    >>> example_seqrecord = SeqRecord(example_seq, id='example')
    >>> new_strip = GenomicStrip(example_seqrecord)
    >>> new_strip.seqrecord.id
    'example'
    >>> str(new_strip.seqrecord.seq)
    'AATTTAATATTGCAGATAGGCAGG'

    The object can be named if desired. If no name is added any reporting or
    graphing method will use the unique ID that is generated automatically
    when the object is created.

    >>> new_strip.add_name('optional name')
    >>> new_strip.name
    'optional name'

    Various types of information can then be attached to the sequence record,
    like lists of segments that are of particular interest.

    >>> segment_1 = Segment(2,5)
    >>> segment_2 = Segment(7,12)
    >>> segment_3 = Segment(13,21)
    >>> segment_list = [segment_1, segment_2, segment_3]
    >>> segment_set = SegmentSet(segment_list, new_strip.seqrecord.id)
    >>> len(new_strip.seg_sets)
    0
    >>> new_strip.add_segment_set(segment_set)
    >>> len(new_strip.seg_sets)
    1
    >>> len(new_strip.seg_sets[0].segments)
    3

    Provenance information can also be added.

    >>>

     """

    def main(self):
        print "This is a GenomicStrip object."

    def __init__(self, seqrecord):
        self.strip_id = create_id('gs')
        self.seqrecord = seqrecord
        self.seg_sets = []
        self.notes = []
        
    def add_name(self, name):
        self.name = name

    def add_segment_set(self, segment_set):
        # store segment data
        self.seg_sets.append(segment_set)
        return self

    def add_note(self, note):
        self.notes.append(note)