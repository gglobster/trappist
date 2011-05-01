__author__ = 'GG'

def td_txt_file_load(filename,tot_items):
    """Load raw info from tab-delimited text file."""
    infile = open(filename, 'r')
    item_count = 0
    while item_count < tot_items :
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

def concatenator(string_list):
    single_string = ""
    for string in string_list :
        single_string = single_string + string
    return single_string

def singleton_list_load(filename,tot_items):
    """Process the first item per line from a tab-delimited text file."""
    rawlines_list = td_txt_file_load(filename,tot_items)
    singletons_list = []
    for rawline in rawlines_list :
        line = strippa_clean(rawline)
        if line[0] == "" :
            break
        singletons_list.append(line[0])
    return singletons_list

def adaptive_list_load(filename,header,columns):
    """Process specified items per line from a tab-delimited text file."""
    rawlines_list = td_txt_file_load(filename,header)
    trimlines_list = []
    for rawline in rawlines_list :
        line = strippa_clean(rawline)
        trim_line = []
        if line[0] == "" :
            break
        for number in columns: # -1 because columns start at 1, not 0
            trim_line.append(line[number-1])
        trimlines_list.append(trim_line)
    return trimlines_list