## script to strip trailing tails from genome file names

import re
from sys import argv
from libs.common import from_dir, ensure_dir
from shutil import copyfile

origin_dir = "data/"+argv[1]+"/"
destin_dir = origin_dir+argv[2]+"/"
file_ext = argv[3]
tail = argv[4]

ensure_dir([destin_dir])

filenames = from_dir(origin_dir, re.compile(r'.*\.'+file_ext))

counter = 0

for filename in filenames:
    # identify strain name
    pattern = re.compile(r'^(.*)'+tail+'\.'+file_ext+'$')
    capture = re.match(pattern, filename)
    # substitute new name
    if capture:
        counter +=1
        new_filename = capture.group(1)+".fas"
        # copy file
        copyfile(origin_dir+filename, destin_dir+new_filename)
        print capture.group(1)

