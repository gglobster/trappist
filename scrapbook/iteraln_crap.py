# starting point and max iterations depend on relative lengths
    if len(seq1) <= len(seq2):
        start = 0
        iter_max = len(seq1)
    else :
        start = len(seq1)-len(seq2)
        iter_max = len(seq2)
    # iterate over the length to find the anchoring position
    iter_count = 0
    scores = []
    while iter_count < iter_max:
        score = 1
        comp_len = 0
        pos1 = start+iter_count
        pos2 = start
        while pos1 < iter_max:
            if seq1[pos1] is seq2[pos2]:
                score +=2
            else:
                score -=1
            pos1 +=1
            pos2 +=1
            comp_len +=1
        # record score
        scores.append(score)
        # test for high score
        if comp_len > 10 and score/comp_len*100 >= 90:
            break
        else:
            iter_count +=1