__author__ = 'GG'

def uniqify(value_list, idfun=None):
    # from http://www.peterbe.com/plog/uniqifiers-benchmark
    """Retain only unique elements while preserving original order."""
    if idfun is None:
        def idfun(x): return x
    seen = {}
    result = []
    for item in value_list:
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

def recursive_bisector(pair_list, size):
    """Bisect coordinate pairs in a list until all are below threshold."""
    # ensure that the values are integers
    max_size = int(size)
    index = 0
    pair_count = len(pair_list)
    while index < pair_count:
        pair = pair_list[index]
        start, end = int(pair[0]), int(pair[1])
        pair_length = end-start
        if pair_length > max_size :
            midpoint = start+(pair_length/2)
            pair_1 = (start, midpoint)
            pair_2 = (midpoint, end)
            pair_list.pop(index)
            pair_list.insert(index, pair_1)
            pair_list.insert(index+1, pair_2)
        else :
            index +=1
        # update length of pair list for loop conditional
        pair_count = len(pair_list)
    return pair_list

def coord_chop(full_length, size, mode):
    """Chop a given length into segments and return coordinate pairs."""
    if size > full_length:
        mode = None
    pair_list = []
    if mode is 'exact_size':
        counter = 0
        while counter < full_length-size:
            new_upper = counter+size
            new_pair = (counter, new_upper)
            pair_list.append(new_pair)
            counter = new_upper
        last_pair = (counter, full_length)
        pair_list.append(last_pair)
    elif mode is 'maxsize_bisect':
        pair_list = [(0, full_length)]
        pair_list = recursive_bisector(pair_list, size)
    elif mode is 'maxsize_divisor':
        divisor = full_length/size +1
        approx_size = full_length/divisor
        counter = 0
        origin = 0
        while counter <= divisor:
            new_upper = origin+approx_size
            new_pair = (origin, new_upper)
            pair_list.append(new_pair)
            origin = new_upper
            counter +=1
        last_pair = (origin, full_length)
        pair_list.append(last_pair)
    elif mode is 'count_divisor':
        seg_count = size # using the size arg to pass the desired seg count
        target_size = full_length/seg_count+1
        pair_list = coord_chop(full_length, target_size, 'exact_size')
    else:
        pair_list = [(0, full_length)]
    return pair_list

def get_region_seq(record,start,stop):
    """Get a single sequence segment by coordinates."""
    segment = record[start:stop]
    segment_seq = segment.seq
    return segment_seq # pure sequence string

def get_seq_subset_by_coords(ori_record, coords_list):
    """Collect a subset of sequence segments by coordinates."""
    from Bio.SeqRecord import SeqRecord
    subset = []
    for coord_pair in coords_list:
        start, stop = coord_pair[0], coord_pair[1]
        segment_seq = get_region_seq(ori_record, start, stop)
        segment_id = ori_record.id+'_'+str(start)+"-"+str(stop)
        new_record = SeqRecord(segment_seq, id=segment_id)
        subset.append(new_record)
    return subset # tuple of sequence records