__author__ = 'GG'

def align_clustal(file_name):
    """Make external call to ClustalW aligner and output a text file."""
    import subprocess
    from Bio.Align.Applications import ClustalwCommandline
    cline = ClustalwCommandline("clustalw", infile=file_name)
    child = subprocess.Popen(str(cline), stdout=subprocess.PIPE, shell=True)
    output, error = child.communicate()
    report = {'output': output, 'error': error}
    # TODO: should set up something to parse ClustalW errors
    return report

def parse_clustal_idstars(file_name):
    """Parse ClustalW output file to measure identity percentage."""
    from analysis.text_manipulation import td_txt_file_load
    import re
    raw_lines = td_txt_file_load(file_name, 3)
    single_line = ''.join(raw_lines)
    idnstar = re.compile(r"\*")
    idnalls = re.findall(idnstar, single_line)
    idntot = len(idnalls)
    return idntot # total number of identical nucleotide positions

def align_muscle(infile_name, outfile_name, log_file):
    """Make external call to Muscle aligner."""
    import subprocess
    from Bio.Align.Applications import MuscleCommandline
    cline = MuscleCommandline(input=infile_name, out=outfile_name, clw=True,
                              loga=log_file, quiet='y') 
    child = subprocess.Popen(str(cline), stdout=subprocess.PIPE, shell=True)
    output, error = child.communicate()
    report = {'output': output, 'error': error}
    # TODO: should set up something to parse MUSCLE errors
    return report

def align_mauve(cline):
    """Make external call to Mauve aligner."""
    import subprocess
    child = subprocess.Popen(str(cline), stdout=subprocess.PIPE, shell=True)
    output, error = child.communicate()
    report = {'output': output, 'error': error}
    # TODO: should set up something to parse Mauve errors
    return report