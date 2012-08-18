from Bio import SeqIO

def write_fasta(filename, seqrecords):
    """Write Fasta file."""
    output_handle = open(filename, 'w')
    count = SeqIO.write(seqrecords, output_handle, 'fasta')
    output_handle.close()
    return count

def write_genbank(filename, seqrecords):
    """Write GenBank file."""
    output_handle = open(filename, 'w')
    counter = SeqIO.write(seqrecords, output_handle, 'genbank')
    output_handle.close()
    return counter