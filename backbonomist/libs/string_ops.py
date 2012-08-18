
def multisplit_finder(g_string, separator):
    """Find multiple split coordinates."""
    origin = 0
    coord = 0
    coord_pairs = []
    while coord < len(g_string):
        coord = g_string.find(separator, origin)
        if coord is -1:
            coord = len(g_string)
        coord_pairs.append((origin, coord))
        origin = coord + len(separator)
    return coord_pairs

# TODO: make sure this correctly handles finding no separators (should return
# TODO: first and last positions of the complete sequence)