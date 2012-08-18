import sys

if sys.argv[1] == "bb_fish":

	from backbonomist import run
	run.main(sys.argv[2:])

elif sys.argv[1] == "bb_map":

	from bb_mapper import run
	run.main(sys.argv[2:])
	

#
#elif sys.argv[1] == "util":
#
#    if sys.argv[2] == "acc_fetch":
#
#        from NGS_utilities import acc_fetcher
#        acc_fetcher.main(sys.argv[1:])
#
#    elif sys.argv[2] == "annot_rem":
#
#        from NGS_utilities import annotate_remote
#        annotate_remote.main(sys.argv[1:])