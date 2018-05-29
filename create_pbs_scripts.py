#!/usr/bin/env python
import sys
import csv, subprocess

parameter_file_full_path = sys.argv[1]
#"/storage/home/fxh140/work/muon_background/veto/PBS/job_params_numu.csv"

with open(parameter_file_full_path, "rb") as csvfile:
    reader = csv.reader(csvfile)
    #('year', 'total_evts', 'num_per_subpart', 'subpart') )
    for job in reader:
        
#        command = """sed -e "s/NUM_DONE_INPUT/0/g" -e "s/RUN_INPUT/{0}/g" -e "s/NUM_PER_SUBPART/{2}/g" -e "s/SUBPART_INPUT/{3}/g" run_stv.sh > /gpfs/group/dfc13/default/dasha/mlarson/pm_trial/nue/pbs/run_stv_nue_{0}_sub{3}.sh""".format(*job)
#        command = """sed -e "s/NUM_DONE_INPUT/0/g" -e "s/RUN_INPUT/{0}/g" -e "s/NUM_PER_SUBPART/{2}/g" -e "s/SUBPART_INPUT/{3}/g" run_stv.sh > /gpfs/group/dfc13/default/dasha/mlarson/pm_trial/numu/pbs/run_stv_numu_{0}_sub{3}.sh""".format(*job)
        command = """sed -e "s/NUM_DONE_INPUT/0/g" -e "s/RUN_INPUT/{0}/g" -e "s/NUM_PER_SUBPART/{2}/g" -e "s/SUBPART_INPUT/{3}/g" run_stv.sh > /gpfs/group/dfc13/default/dasha/mlarson/pm_trial/nutau/pbs/run_stv_nutau_{0}_sub{3}.sh""".format(*job)
#        command = """sed -e "s/NUM_DONE_INPUT/0/g" -e "s/RUN_INPUT/{0}/g" -e "s/NUM_PER_SUBPART/{2}/g" -e "s/SUBPART_INPUT/{3}/g" run_stv.sh > /gpfs/group/dfc13/default/dasha/mlarson/pm_trial/corsika/pbs/run_stv_corsika_{0}_sub{3}.sh""".format(*job)
        
        print command # Uncomment this line when testing to view the qsub command

        # Comment the following 3 lines when testing to prevent jobs from being submitted
        exit_status = subprocess.call(command, shell=True)
        if exit_status is 1:  # Check to make sure the job submitted
            print "script run_estes_muongun_{0}_sub{3}.sh failed to generate".format(*job)
print "Done!"
