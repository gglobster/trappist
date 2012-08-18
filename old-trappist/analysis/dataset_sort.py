__author__ = 'GG'

def reorder_matches_by_contig(q_matches):
    """Reorder query match sets by contig."""
    cc_matches = {}
    for q_match in q_matches:
        for match in q_match['matches']:
            c_match = {'query_id': q_match['query_id'],
                       'details': match['details']}
            if match['contig_id'] in cc_matches.keys():
                cc_matches[match['contig_id']].append(c_match)
            else:
                cc_matches[match['contig_id']] = [c_match]
    return cc_matches