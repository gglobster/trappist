__author__ = 'GG'


### Sequence feature operations ###

def feat_collect(infile, feat_mode):
    """Collect a subset of features from a genbank file."""
    from analysis.seqfile_ops import load_genbank
    gb_record = load_genbank(infile)
    feat_list = gb_record.features
    collected = []
    # establish collection parameters
    types_list = feat_mode['types'] # default entry is ('CDS')
    tags_dict = feat_mode['tags'] # default is an empty dictionary
    # start collecting features
    for feature in feat_list:
        if feature.type in types_list:
            if len(tags_dict.keys()) is 0:
                collected.append(feature)
            else:
                for tag_key in tags_dict.keys():
                    if tag_key in feature.qualifiers:
                        feat_value = feature.qualifiers.get(tag_key)
                        if feat_value[0] in tags_dict[tag_key]:
                            collected.append(feature)
                        else: pass
                    else: pass
        else: pass
    ## consider adding some info to the log
    return collected

def feature_coords(features):
    """Extract coordinates from a list of sequence features."""
    coords_list = []
    for feature in features:
        coord_start = feature.location.nofuzzy_start
        coord_end = feature.location.nofuzzy_end
        coord_pair = (coord_start, coord_end)
        coords_list.append(coord_pair)
    ## consider adding some info to the log
    return coords_list


### Sequence record operations ###

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
        start, stop = int(coord_pair[0]), int(coord_pair[1])
        segment_seq = get_region_seq(ori_record, start, stop)
        segment_id = ori_record.id+'_'+str(start)+"-"+str(stop)
        new_record = SeqRecord(segment_seq, id=segment_id)
        subset.append(new_record)
    return subset # tuple of sequence records


### Coordinates chop operations ###

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
