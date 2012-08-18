# script to generate a genome set file for bb_mapper from dir contents

import re
from sys import argv
from libs.common import from_dir, load_fasta, load_multifasta, load_genbank

data_dir = "data/"+argv[1]+"/"
seq_dir = data_dir+argv[2]+"/"
py_out = data_dir+argv[3]+"_set.py"
min_size = argv[4]

set_lines = ["all = ["]

filenames = from_dir(seq_dir, re.compile(r'.*\..*'))

counter = 1

for filename in filenames:

    print filename,

    while True:

        if filename.find(".gbk") > 0:
            # process genbank
            try:
                record = load_genbank(seq_dir+filename)
            except IOError:
                print "failed to load Genbank file"
                break
            except Exception:
                print "failed to handle Genbank file"
                break
            else:
                print "...",
                seq_format = 'gbk'

        elif filename.find(".fas") > 0:
            # process fasta (for mfas, load first record)
            try:
                record = load_fasta(seq_dir+filename)
            except IOError:
                print "failed to load Fasta file as single-record file"
                break
            except Exception:
                try:
                    record = load_multifasta(seq_dir+filename)[0]
                except IOError:
                    print "failed to load Fasta file as multi-record file"
                    break
                except Exception:
                    print "failed to handle Fasta file"
                    break
            print "...",
            seq_format = 'fas'

        else:
            # reject as bad format
            print "invalid file format"
            break

        if len(record) < int(min_size):
            print "too short --", len(record)
            break

        line = "".join(["\t{'name': '", record.id,
                        "', 'file': '", filename,
                        "', 'input': '", seq_format,
                        "', 'order': ", str(counter),
                        ", 'nudge': 0, 'offset': 0,'ignore': (0, 0)},"])
                        # TODO: edit to fit backbonomist or bbmapper
        set_lines.append(line)
        counter +=1

        print "OK"
        break

set_lines.append("]")

open(py_out, 'w').write("\n".join(set_lines))


