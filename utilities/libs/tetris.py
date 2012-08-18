
def segment_finder(reference, vectorDB, keyword, namebase):
    """Find conserved segments in a set of vectorized sequences."""

    segment_DB = {}
    segment_id = None
    seg_init = True

    for symbol in reference:

        synteny = False
        ref_pos = reference.index(symbol)

        while True:

            for v_key in vectorDB:

                synteny = False
                query = vectorDB[v_key][keyword]

                try:
                    position1 = query.index(reference[ref_pos])
                except IndexError:
                    break
                except ValueError:
                    break
                else:
                    if position1:
                        try:
                            position2 = query.index(reference[ref_pos+1])
                        except IndexError:
                            break
                        except ValueError:
                            break
                        else:
                            if position2:
                                if position2 == position1+1:
                                    synteny = True
                            else:
                                break
                    else:
                        break

            ref_pos +=1
            break

        if synteny:
            if seg_init:
                segment_id = namebase+str(len(segment_DB)+1)
                # not sure what's up with indexing but this works
                segment_DB[segment_id] = [reference[ref_pos-1],
                                          reference[ref_pos]]
                # switch off initialization
                seg_init = False
            else:
                segment_DB[segment_id].append(reference[ref_pos])
        else:
            seg_init = True
    return segment_DB