import csv
import sys
import event_counter
import os

no_evt = sys.argv[2]
f = open(sys.argv[1], 'wt')
try:
    writer = csv.writer(f,quoting=csv.QUOTE_ALL)
    for run in range(0,1000):
        run=str(run)
        run_pad=run.zfill(6)
        print "run_pad = ", run_pad
        input_file = "/gpfs/group/dfc13/default/dasha/mlarson/L2/corsika_all/Level2_IC86.2012_corsika.011808.%s.i3.bz2"% run_pad
#        input_file = "/gpfs/group/dfc13/default/dasha/mlarson/L2/nue/genie_ic.12640.%s.i3.gz"  % run_pad
#        input_file = "/gpfs/group/dfc13/default/dasha/mlarson/L2/numu/genie_ic.14640.%s.i3.gz"  % run_pad
#        input_file = "/gpfs/group/dfc13/default/dasha/mlarson/L2/nutau/genie_ic.16640.%s.i3.gz"  % run_pad
        if os.path.isfile(input_file):
            print "input_file = ", input_file
            total_num_evts = event_counter.main(input_file)
            print total_num_evts
            n, mod =divmod(total_num_evts, int(no_evt))
            print "n, mod  = ", n, mod 
            for i in range(0, n+1):
                writer.writerow( (run_pad, total_num_evts, no_evt, i ) )
finally:
    f.close()

#print open(sys.argv[1], 'rt').read()
