__author__ = 'GG'

def td_txt_file_load(filename, skip_items):
    """Load raw info from tab-delimited text file."""
    infile = open(filename, 'r')
    item_count = 0
    # option to skip a few lines at the top
    while item_count < skip_items :
        infile.readline()
        item_count += 1
    rawlines_list = infile.readlines()
    infile.close()
    return rawlines_list

def strippa_clean(rawline):
    """Clean up raw data line from td file read-in into a line tuple."""
    stripline = rawline.strip("\n")
    stripline = stripline.strip("\r")
    line_tuple = stripline.split("\t")
    return line_tuple

def adaptive_list_load(filename, skip_items, columns):
    """Process specified items per line from a tab-delimited text file."""
    rawlines_list = td_txt_file_load(filename, skip_items)
    trimlines_list = []
    for rawline in rawlines_list :
        line = strippa_clean(rawline)
        trim_line = []
        if line[0] == "" :
            break
        for number in columns:
            trim_line.append(line[number])
        trimlines_list.append(trim_line)
    return trimlines_list