__author__ = 'GG'

# to switch between gbk and fasta, use
# count = SeqIO.convert(infile, "genbank", outfile, "fasta")
# and vice-versa

# should add a function to manage directory handling

def seqfile_format(filename):
    """Identify file format based on file extension."""
    # must be a better way to do this
    import re
    gb_ext = re.compile(r"(.*)\.gb")
    fa_ext = re.compile(r"(.*)\.fa")
    sq_ext = re.compile(r"(.*)\.seq")
    genbank = gb_ext.search(filename)
    fasta = fa_ext.search(filename)
    seq = sq_ext.search(filename)
    if genbank :
        format = 'genbank'
        core_gb = gb_ext.match(filename)
        name = core_gb.group(1)
    elif fasta :
        format = 'fasta'
        core_fa = fa_ext.match(filename)
        name = core_fa.group(1)
    elif seq :
        format = 'seq'
        core_sq = sq_ext.match(filename)
        name = core_sq.group(1)
    else :
        format = 'unidentified'
        name = filename
    return format, name

def ensure_fasta(filename):
    """Check file format; if it's GenBank, convert to Fasta."""
    # deprecate in favor of surefmt_load
    from Bio import SeqIO
    format, name = seqfile_format(filename)
    if format == 'genbank':
        outfile = name+".fas"
        SeqIO.convert(filename, 'genbank', outfile, 'fasta')
    else:
        outfile = filename
    return outfile

def surefmt_load(filename, need_format, alphabet):
    """Load file with format check; convert if needed."""
    from Bio import SeqIO
    seq_record, rec_type = load_agnostic(filename)
    if rec_type is not need_format:
        if need_format is 'genbank': ext = ".gbk"
        elif need_format is 'fasta': ext = ".fas"
        else: ext = ".txt"
        outfile = seq_record.id+ext
        try: SeqIO.convert(filename, rec_type, outfile, need_format, alphabet)
        except: raise
        else: seq_record, rec_type = load_agnostic(outfile)
    else:
        pass
    return seq_record

def write_fasta(filename, seqrecords):
    """Write Fasta file."""
    from Bio import SeqIO
    output_handle = open(filename, 'w')
    count = SeqIO.write(seqrecords, output_handle, 'fasta')
    output_handle.close()
    return count

def write_genbank(filename, seqrecords):
    """Write GenBank file."""
    from Bio import SeqIO
    output_handle = open(filename, 'w')
    counter = SeqIO.write(seqrecords, output_handle, 'genbank')
    output_handle.close()
    return counter

def load_agnostic(seqfile):
    """Load single-record file of unspecified format."""
    # present version is single-record fasta or genbank only
    # TODO: expand options
    while True:
        try: seq_record = load_fasta(seqfile)
        except: print "not single-record fasta"
        else: 
            rec_type = 'fasta'
            print "found a fasta file"
            return seq_record, rec_type
            break
        try: seq_record = load_genbank(seqfile)
        except: print "not genbank"
        else: 
            rec_type = 'genbank'
            print "found a genbank file"
            return seq_record, rec_type
            break
        raise Exception("Cannot open query file!")

def load_fasta(seqfile):
    """Load single-record Fasta file."""
    from Bio import SeqIO
    input_handle = open(seqfile, 'rU')
    fasta_record = SeqIO.read(input_handle, 'fasta')
    assert fasta_record.id 
    input_handle.close()
    return fasta_record

def load_multifasta(seqfile):
    """Load multiple records from Fasta file."""
    from Bio import SeqIO
    input_handle = open(seqfile, 'rU')
    multifasta = SeqIO.parse(input_handle, 'fasta')
    fasta_list = list(multifasta)
    for record in fasta_list:
        assert record.id
    return fasta_list

def load_genbank(seqfile):
    """Load single-record GenBank file."""
    from Bio import GenBank
    parser = GenBank.FeatureParser()
    input_handle = open(seqfile, 'rU')
    gb_record = parser.parse(input_handle)
    assert gb_record.id 
    input_handle.close()
    return gb_record

def compile_multifasta_from_folders(base_names_list,data_dir,folder_ext,
                                    file_ext,number_base):
    """Compile fasta records from folders into single multifasta files."""
    # TODO: rework and make unit test
    from os import path
    glob_count  = 0
    for base_name in base_names_list:
        print "processing", base_name
        #print  base_name, base_name+".fas", "makeDB"
        dir_path = data_dir + base_name + folder_ext +"/"
        counter = number_base
        fas_records = []
        while 1 :
            seq_infile = dir_path + base_name + str(counter) + file_ext
            if not path.exists(seq_infile):
                print "no such file"
                break
            else:
                record = load_fasta(seq_infile)
                print record.id
                fas_records.append(record)
                counter += 1
        seq_outfile = data_dir + base_name +".fas"
        fas_count = write_fasta(seq_outfile,fas_records)
        glob_count = glob_count + fas_count
    return glob_count
