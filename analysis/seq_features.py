__author__ = 'GG'

def feat_collect(infile, feat_mode):
    """Collect a subset of features from a genbank file."""
    from analysis.sequence_file_ops import load_genbank
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
                        if feat_value in tags_dict[tag_key]:
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
