from loaders import load_fasta, load_genbank
from writers import write_genbank, write_fasta
from Bio.Alphabet import generic_dna

def fas2gbk(fas_file):
    """Convert a FastA file to Genbank format."""
    record = load_fasta(fas_file)
    gbk_file = fas_file[:fas_file.find('.fas')]+'.gbk'
#    record.name = rec_name
#    record.id = rec_name
    record.seq.alphabet = generic_dna
    write_genbank(gbk_file, record)
    return gbk_file

def gbk2fas(gbk_file):
    """Convert a Genban file to kFastA format."""
    record = load_genbank(gbk_file)
    fas_file = gbk_file[:gbk_file.find('.gbk')]+'.fas'
#    record.name = rec_name
#    record.id = rec_name
    write_fasta(fas_file, record)
    return fas_file