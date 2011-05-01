__author__ = 'GG'

def align_clustal(file_name):
    """Make external call to ClustalW aligner and output a text file."""
    import sys,subprocess
    from Bio.Align.Applications import ClustalwCommandline
    cline = ClustalwCommandline("clustalw", infile=file_name)
    child = subprocess.Popen(str(cline), stdout=subprocess.PIPE,
                             shell=(sys.platform!="win32"))
    output, error = child.communicate()
    # should do something to parse ClustalW errors
    return output

def parse_clustal_idstars(file_name):
    """Parse ClustalW output file to measure identity percentage."""
    from analysis.text_manipulation import td_txt_file_load
    import re
    raw_lines = td_txt_file_load(file_name,3)
    single_line = ''.join(raw_lines)
    idnstar = re.compile(r"\*")
    idnalls = re.findall(idnstar,single_line)
    idntot = len(idnalls)
    return idntot # total number of identical nucleotide positions

def align_muscle(infile_name,outfile_name):
    """Make external call to Muscle aligner."""
    import sys,subprocess
    from Bio.Align.Applications import MuscleCommandline
    cline = MuscleCommandline(input=infile_name, out=outfile_name,
                              clw=True, loga='muscle_log.txt',
                              verbose='yes') #, quiet='no'
                              # enable quiet/verbose switch at config time
    child = subprocess.Popen(str(cline), stdout=subprocess.PIPE,
                             shell=(sys.platform!="win32"))
    output, error = child.communicate()
    # should do something to parse MUSCLE errors
    return output

def align_mauve(cline):
    """Make external call to Mauve aligner."""
    import sys,subprocess
    child = subprocess.Popen(str(cline), stdout=subprocess.PIPE,
                             shell=(sys.platform!="win32"))
    output, error = child.communicate()
    # should do something to parse Mauve errors
    return output