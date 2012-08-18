def FastqJointIterator(handleA, handleB):
    """Iterate over Fastq records from two separate files as string tuples.

    Modified from Bio.SeqIO.QualityIO import FastqGeneralIterator to read
    from two separate files in parallel. Useful for shuffle-combining
    Illumina read sets.

    """
    handleA_readline = handleA.readline
    handleB_readline = handleB.readline
    while True:
        lineA = handleA_readline()
        lineB = handleB_readline()
        if lineA == "" or lineB == "": return
        if lineA[0] == "@" or lineB[0] == "@":
            break
    while True:
        title_lines = []
        seq_strings = []
        quality_strings = []
        count = 0
        while count < 2 :
            if lineA[0] != "@" and lineB[0] != "@":
                raise ValueError("Bad formatting of record start lines.")
            title_lineA = lineA[1:].rstrip()
            title_lineB = lineB[1:].rstrip()
            title_lines.append((title_lineA, title_lineB))
            seq_stringA = handleA_readline().rstrip()
            seq_stringB = handleB_readline().rstrip()
            seq_strings.append((seq_stringA, seq_stringB))
            while True:
                lineA = handleA_readline()
                lineB = handleB_readline()
                if not lineA:
                    raise ValueError("End of file without quality info (A).")
                if not lineB:
                    raise ValueError("End of file without quality info (B).")
                if lineA[0] == "+" and lineB[0] == "+":
                    second_titleA = lineA[1:].rstrip()
                    second_titleB = lineB[1:].rstrip()
                    if second_titleA and second_titleA != title_lineA:
                        raise ValueError("Seq and qual captions differ (A).")
                    if second_titleB and second_titleB != title_lineB:
                        raise ValueError("Seq and qual captions differ (B).")
                    break
                seq_stringA += lineA.rstrip() #removes trailing newlines
                seq_stringB += lineB.rstrip() #removes trailing newlines
            if " " in seq_stringA or "\t" in seq_stringA:
                raise ValueError("Whitespace not allowed in sequence (A).")
            if " " in seq_stringB or "\t" in seq_stringB:
                raise ValueError("Whitespace not allowed in sequence (B).")
            seq_lenA = len(seq_stringA)
            seq_lenB = len(seq_stringB)
            quality_stringA = handleA_readline().rstrip()
            quality_stringB = handleB_readline().rstrip()
            quality_strings.append((quality_stringA, quality_stringB))
            while True:
                lineA = handleA_readline()
                lineB = handleB_readline()
                if not lineA or not lineB: break #end of file
                if lineA[0] == "@":
                    if len(quality_stringA) >= seq_lenA :
                        break
                quality_stringA += lineA.rstrip()
                if lineB[0] == "@":
                    if len(quality_stringB) >= seq_lenB:
                        break
                quality_stringB += lineB.rstrip()
            if seq_lenA != len(quality_stringA):
                raise ValueError("Lengths of seq and qual values differs "
                                 " for %s (%i and %i) (A)." \
                                 % (title_lineA, seq_lenA,
                                    len(quality_stringA)))
            if seq_lenB != len(quality_stringB):
                raise ValueError("Lengths of seq and qual values differs "
                                 " for %s (%i and %i) (B)." \
                                 % (title_lineB, seq_lenB,
                                    len(quality_stringB)))
            count +=1
        yield (title_lines, seq_strings, quality_strings)
        if not lineA or not lineB: return #StopIteration at end of file
    assert False, "Should not reach this line"

def combine_illumina(dataset): # without demuxing
    """Combine forward and reverse read sets from Illumina.

    Consolidate both read sets in a single FastQ file, with each forward read
    immediately followed by the corresponding reverse read. This is meant to
    replace the "shuffle" Perl script provided with the Velvet assembler.

    """
    print " ", dataset['run_id']
    # identify inputs and outputs
    ori_root = dirs['ori_data_dir']
    run_nick = dataset['run_id']
    fwd_file = ori_root+run_nick+"/"+dataset['source_fwd']
    rev_file = ori_root+run_nick+"/"+dataset['source_rev']
    combined_file = dirs['combined_dir']+run_nick+".txt"
    # start shuffling
    co_out = open(combined_file, 'w')
    pair_count = 0
    for titles, seqs, quals in FastqJointIterator(open(fwd_file), 
                                                  open(rev_file)) :
        F_title, R_title = titles
        F_seq, R_seq = seqs
        F_qual, R_qual = quals
        # output to combined file
        co_out.write("@%s\n%s\n+\n%s\n" % (F_title, F_seq, F_qual))
        co_out.write("@%s\n%s\n+\n%s\n" % (R_title, R_seq, R_qual))
        # increment counter
        pair_count +=1
        # report on the progress
        if pair_count%100000==0:
            print "\t", pair_count, "read pairs processed"
    co_out.close()