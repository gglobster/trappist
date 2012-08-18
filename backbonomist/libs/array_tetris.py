import numpy as np

def extract_nonzero(array_o, array_d):
    """Extract rows that are non-zero into an array of same shape."""
    for row in array_o:
        a, b, c, d = row
        if a == 0 or c == 0:
            pass
        else:
            array_d = np.append(array_d, row)
    return array_d

def clump_rows(array, dist, mtype):
    """Collapse rows of coordinates that are close and in same orientation."""
    # sort array to ensure optimal clumping
    array = np.sort(array, order=mtype[0][0])
    # perform clumping
    r_ix = 1
    new_len = len(array)
    while r_ix < new_len:
        xa1, xb1, xc1, xd1 = array[r_ix-1]
        xa2, xb2, xc2, xd2 = array[r_ix]
        # clump if test conditions are fulfilled
        if float(xa1)/float(xc1) > 0 : # same sign
            gap1 = abs(abs(xa2)-abs(xb1))
            gap2 = abs(abs(xc2)-abs(xd1))
            if gap1 < dist and gap2 < dist and 0.75 < gap1/(gap2+1) < 1.3:
                xa,xb = xa1,xb2
                xc,xd = xc1,xd2
                new_row = xa, xb, xc, xd
                # delete the original rows
                array = np.delete(array, (r_ix-1, r_ix), 0)
                # insert the updated row
                array = np.insert(array, r_ix-1, new_row)
            else :
                r_ix += 1
        elif float(xa1)/float(xc1) < 0 : # different sign
            gap1 = abs(abs(xa1)-abs(xb2))
            gap2 = abs(abs(xc2)-abs(xd1))
            if gap1 < dist and gap2 < dist and 0.75 < gap1/(gap2+1) < 1.3:
                xa,xb = xa2,xb1
                xc,xd = xc1,xd2
                new_row = xa,xb,xc,xd
                # delete the original rows
                array = np.delete(array, (r_ix-1, r_ix), 0)
                # insert the updated row
                array = np.insert(array, r_ix-1, new_row)
            else :
                r_ix += 1
        else : # this should not happen
            print "WTF? This shouldn't happen."
            r_ix += 1
        new_len = len(array)
    return array

def chop_rows(uncut_array, max_size, chop_mode, mtype):
    """Chop rows of coordinates that span too large segments."""
    # set up an array stub to receive new rows
    all_rows_list = []
    for xa, xb, xc, xd in uncut_array:
        size1 = abs(abs(xa)-abs(xb))
        pair_list1 = coord_chop(size1, max_size, chop_mode)
        size2 = abs(abs(xc)-abs(xd))
        pair_list2 = coord_chop(size2, max_size, chop_mode)
        try:
            assert len(pair_list1) == len(pair_list2)
        except AssertionError:
            # edge case, cancel chopping for this segment
            row = xa, xb, xc, xd
            all_rows_list.append(row)
        else:
            # denormalize coordinates
            dn_pairs1 = []
            for a, b in pair_list1:
                if xa < 0: pair = (xa-a, xa-b)
                else: pair = (xa+a, xa+b)
                dn_pairs1.append(pair)
            dn_pairs2 = []
            for c, d in pair_list2:
                if xc < 0: pair = (xc-c, xc-d)
                else: pair = (xc+c, xc+d)
                dn_pairs2.append(pair)
            # consolidate to row format
            count = 0
            while count < len(dn_pairs1):
                row = (dn_pairs1[count][0], dn_pairs1[count][1],
                       dn_pairs2[count][0], dn_pairs2[count][1])
                all_rows_list.append(row)
                count +=1
    # make a new array from the list of rows
    chop_array = np.array(all_rows_list, dtype=mtype)
    chop_array = np.sort(chop_array, order=mtype[0][0])
    return chop_array

def coord_chop(length, max_size, mode):
    """Chop a given length into segments and return coordinate pairs."""
    if length < max_size:
        mode = None
    pair_list = []
    if mode is 'exact_size':
        counter = 0
        while counter < length-max_size:
            new_upper = counter+max_size
            new_pair = (counter, new_upper)
            pair_list.append(new_pair)
            counter = new_upper
        last_pair = (counter, length)
        pair_list.append(last_pair)
    elif mode is 'maxsize_bisect':
        pair_list = [(0, length)]
        pair_list = recursive_bisector(pair_list, max_size)
    elif mode is 'maxsize_divisor':
        divisor = length/max_size +1
        approx_size = length/divisor
        counter = 0
        origin = 0
        while counter <= divisor:
            new_upper = origin+approx_size
            new_pair = (origin, new_upper)
            pair_list.append(new_pair)
            origin = new_upper
            counter +=1
        last_pair = (origin, length)
        pair_list.append(last_pair)
    elif mode is 'count_divisor':
        seg_count = max_size # using the size arg to pass the desired seg count
        target_size = length/seg_count+1
        pair_list = coord_chop(length, target_size, 'exact_size')
    else:
        pair_list = [(0, length)]
    return pair_list

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

def get_anchor_loc(quad_array):
    """Determine which segment is largest and where it hits."""
    # set up an array stub to receive size and location data
    sl_array = np.ones(1, dtype=[('size', 'i4'),
                                 ('start', 'i4'),
                                 ('end', 'i4'),
                                 ('orient', 'i2')])
    # collect data
    for row in quad_array:
        a, b, c, d = row
        size = -abs(abs(a)-abs(b)) # size negative'd so biggest shows up first
        loc_start = abs(a)
        loc_end = abs(b)
        orient = abs(a)/a
        data = size, loc_start, loc_end, orient
        sl_array = np.insert(sl_array, 0, data)
    # evaluate which segment is largest
    sl_array = np.sort(sl_array, order='size')
    anchor = sl_array[0]
    return anchor

def offset_q2r_coords(data_array, q_len, value_pair, segtype):
    """Offset segment coordinates in array to specified values.

    This is tailored for pairwise alignments in which the reference and
    query sequences are not treated equally. The reference is treated as a
    linear molecule and its coordinates (first pair in array quintet) get
    nudged to the right. The query is treated as a circular molecule and its
    coordinates (second pair in array quintet) are looped around the zero
    position.
    """
    # set up an array stub to receive new rows
    all_rows_list = []
    # run offsetting
    r_off, q_off = value_pair
    for xa, xb, xc, xd, idp in data_array:
        xao = nudge_coord(xa, r_off)
        xbo = nudge_coord(xb, r_off)
        xco = offset_coord(xc, q_len, q_off)
        xdo = offset_coord(xd, q_len, q_off)
        offset_row = (xao, xbo, xco, xdo, idp)
        all_rows_list.append(offset_row)
    # make a new array from the list of rows
    offset_array = np.array(all_rows_list, dtype=segtype)
    return offset_array

def offset_coord(coord, length, offset) :
    """Calculate coordinate offset."""
    # calculate linear coordinate change
    coft = (abs(coord)-abs(offset)+1)
    if coft < 0 :
        coff = length + coft
    else :
        coff = coft
    # invert sequence if offset is negative
    if offset < 0 :
        cofi = length - coff
    else :
        cofi = coff
    # recall original sign
    if coord < 0 :
        sign0 = -1
    else :
        sign0 = 1
    # set final sign
    coff7 = sign0 * cofi
    return coff7

def nudge_coord(coord, offset) :
    """Calculate coordinate offset."""
    # calculate linear coordinate change
    coft = (abs(coord)+abs(offset)+1)
    # recall original sign
    if coord < 0 :
        sign0 = -1
    else :
        sign0 = 1
    # set final sign
    coff7 = sign0 * coft
    return coff7

def shade_split(xa, xb, xc, xd, q_len):
    """Split shaded areas that sit across the map origin.

    This assumes that xa -> xb is always non-split since it is on the
    reference, which is treated as linear, not circular. It also assumes that
    only reference-side segments can be of negative sign.
    """
    if xa > 0:
        xa1, xb1, xc1, xd1 = xa, xa+q_len-xc, xc, q_len
        xa2, xb2, xc2, xd2 = xa+q_len-xc, xb, 1, xd
    else: # xa <0
        xa1, xb1, xc1, xd1 = xa, xa-xd, 1, xd
        xa2, xb2, xc2, xd2 = xa-xd, xb, xc, q_len
    coords1 = xa1, xb1, xc1, xd1
    coords2 = xa2, xb2, xc2, xd2
    return coords1, coords2

def coord_flipper(axr,bxr) :
    """Flip negative coordinates."""
    if axr < 0 :
        ax,bx = abs(bxr),abs(axr)
    else :
        ax,bx = axr,bxr
    return ax,bx




























