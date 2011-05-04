__author__ = 'GG'

def uniqify(seq,idfun=None):
    # from http://www.peterbe.com/plog/uniqifiers-benchmark
    """Retain only unique elements while preserving original order."""
    if idfun is None:
        def idfun(x): return x
    seen = {}
    result = []
    for item in seq:
        marker = idfun(item[0])
        if marker in seen: continue
        seen[marker] = 1
        result.append(item)
    return result

def ensure_dir(dir_path):
    """Check that the directory exists; if not, create it."""
    from os import path, makedirs
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
    return abs_path, report

def coord_chop(length, size, mode):
    """Chop a given length into segments and return coordinate pairs."""
    if mode is 'exact_size':
        pair_list = []
        counter = 0
        while counter < length-size:
            new_upper = counter+size
            new_pair = (counter, new_upper)
            pair_list.append(new_pair)
            counter = new_upper
        last_pair = (counter, length)
        pair_list.append(last_pair)
        print pair_list # INFO logging
    elif mode is 'max_size':
        # chop recursively in two like in Bacon
        pair_list = (0,length)  ## TBD
    else:
        pair_list = (0,length)
    return pair_list

def get_region_seq(record,start,stop):
    """Get sequence segment."""
    segment = record[start:stop]
    segment_seq = segment.seq
    return segment_seq # pure sequence string

def get_seq_subset_by_coords(ori_record, coords_list):
    """Collect a subset of sequence segments based on coordinates."""
    from Bio.SeqRecord import SeqRecord
    subset = []
    for coord_pair in coords_list:
        start, stop = coord_pair[0], coord_pair[1]
        segment_seq = get_region_seq(ori_record, start, stop)
        segment_id = str(start)+"_"+str(stop)
        new_record = SeqRecord(segment_seq, id=segment_id)
        subset.append(new_record)
    return subset # tuple of sequence records