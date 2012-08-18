from Bio import SeqIO
from os import path, makedirs

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

def ensure_dir(dir_list):
    """Check that the directory exists; if not, create it."""
    for dir_path in dir_list:
        abs_path = path.abspath(dir_path)
        if not path.exists(abs_path):
            try: makedirs(abs_path)
            except Exception as message:
                status = 1
                # TODO: make graceful fail or request input if interactive mode
            else:
                message = 'created path'
                status = 0
        else:
            message = 'path exists'
            status = 0
        report = {'message': message, 'status': status} 

