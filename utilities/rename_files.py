## script to rename and copy sets of files

import re
from sys import argv
from libs.common import from_dir, ensure_dir
from shutil import copyfile

origin_dir = "data/"+argv[1]
destin_dir = "data/"+argv[2]+"/"
prefix = argv[3]
postfix = argv[4]
sub_base = argv[5]

ensure_dir([destin_dir])

filenames = from_dir(origin_dir, re.compile(r'.*\.fas.*'))

counter = 0

for filename in filenames:
    # identify strain name
    pattern = re.compile(r'^'+prefix+'(.*)'+postfix+'$')
    capture = re.match(pattern, filename)
    # substitute new name
    if capture:
        counter +=1
        new_filename = sub_base+"_"+str(counter)+".fas"
        # copy file
        copyfile(origin_dir+"/"+filename, destin_dir+new_filename)
        print capture.group(1), sub_base+"_"+str(counter)

