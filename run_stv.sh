#!/bin/sh
#PBS -W group_list=dfc13_collab
#PBS -l nodes=1:ppn=1
#PBS -l mem=4gb
#PBS -l walltime=8:00:00
#PBS -j oe
#PBS -o /storage/home/dup193/work/estes/output/test.log
#PBS -N nutau
#PBS -m n

start=`date +%s`

#==============================================================================
# Setup environment similar to what's done in .bashrc
#
# (Also clears out PATH and PYTHONPATH so that user's settings don't interfere
# with this script)
#==============================================================================

UHOME=~dup193
export PARROT_RUN="$UHOME/cctools/bin/parrot_run"
export PARROT_CVMFS_ALIEN_CACHE="/gpfs/group/dfc13/default/cache/"
export CVMFS="/cvmfs/icecube.opensciencegrid.org"
ENV_SHELL="/gpfs/group/dfc13/default/dasha/icerec_V05-01-05/build/env-shell.sh"
SCRIPT="/storage/home/dup193/work/estes/new_estes/STV/STV.py"
#SCRIPT="/storage/home/fxh140/work/muon_background/veto/ESTES_run.py"

#==============================================================================
# Define functions whereby CVMFS is accessible
# via the Parrot user-space tool (and is intended--but not tested--to also work
# with a proper native installation of CVMFS)
#==============================================================================

function init_i3_env () {
	py2_vx=$1
	shift
	env_shell=$1
	shift
	if [ -x "$CVMFS/${py2_vx}/setup.sh" ]
	then
		eval $( $CVMFS/${py2_vx}/setup.sh )
	else
		echo "ERROR: $CVMFS/${py2_vx}/setup.sh not found or not executable!"
	fi

	if [ -x "${env_shell}" ]
	then
		${env_shell} "$*"
	else
		echo "ERROR: ${env_shell} not found or not executable!"
	fi
} 
export -f init_i3_env

function i3_run () {
	if [ ! -d "$CVMFS" ]
	then
		if [ -x "`which cvmfs_config 2>/dev/null`" ]
		then
			cvmfs_config probe
			/bin/sh --norc -c "'init_i3_env' $*"
		elif [ ! -d "$CVMFS" -a -x "$PARROT_RUN" ]
		then
			HTTP_PROXY="cache01.hep.wisc.edu:3128" "$PARROT_RUN" /bin/sh --norc -c "'init_i3_env' $*"
		fi
	fi
}
export -f i3_run



#outdir=/gpfs/group/dfc13/default/dasha/mlarson/output_1/nue/
#outdir=/gpfs/group/dfc13/default/dasha/mlarson/output_1/numu/
outdir=/gpfs/group/dfc13/default/dasha/mlarson/output_1/nutau/
#outdir=/gpfs/group/dfc13/default/dasha/mlarson/output_1/corsika/

if [ ! -d "$DIRECTORY" ]; then
    mkdir -p ${outdir}
fi
cd ${outdir}

#num_done=0  # no. of events already done, start running from the num_done+1 th event, e.g. num_done =2, then scripts only run from the 3th event                                                          
#num_events=100        # no. of events per subpart                                                   
#subpart=0                                                                                           


num_done=NUM_DONE_INPUT                                                                               
run_pad=RUN_INPUT                                                                                     
subpart=SUBPART_INPUT
echo 'run: '$run_pad
echo 'subpart: '$subpart
num_events=NUM_PER_SUBPART  
echo "num_events per subpart: "$num_events
num_skip=$(($subpart*$num_events + $num_done))
#num_skip=`expr $subpart \* $num_events \+ $num_done`
echo 'skip: '$num_skip

#i3_run py2-v2 "$ENV_SHELL" 'python' "$SCRIPT" -i /gpfs/group/dfc13/default/dasha/mlarson/L2/GeoCalibDetectorStatus_2013.56429_V1_Modified.i3.gz /gpfs/group/dfc13/default/dasha/mlarson/L2/nue/genie_ic.12640.$run_pad.i3.gz -o /gpfs/group/dfc13/default/dasha/mlarson/output_1/nue/Pm_genie_ic.12640.$run_pad.part$subpart --sk $num_skip --ne $num_events 

#i3_run py2-v2 "$ENV_SHELL" 'python' "$SCRIPT" -i /gpfs/group/dfc13/default/dasha/mlarson/L2/GeoCalibDetectorStatus_2013.56429_V1_Modified.i3.gz /gpfs/group/dfc13/default/dasha/mlarson/L2/numu/genie_ic.14640.$run_pad.i3.gz -o /gpfs/group/dfc13/default/dasha/mlarson/output_1/numu/Pm_genie_ic.14640.$run_pad.part$subpart --sk $num_skip --ne $num_events 

i3_run py2-v2 "$ENV_SHELL" 'python' "$SCRIPT" -i /gpfs/group/dfc13/default/dasha/mlarson/L2/GeoCalibDetectorStatus_2013.56429_V1_Modified.i3.gz /gpfs/group/dfc13/default/dasha/mlarson/L2/nutau/genie_ic.16640.$run_pad.i3.gz -o /gpfs/group/dfc13/default/dasha/mlarson/output_1/nutau/Pm_genie_ic.16640.$run_pad.part$subpart --sk $num_skip --ne $num_events 

#i3_run py2-v2 "$ENV_SHELL" 'python' "$SCRIPT" -i /gpfs/group/dfc13/default/dasha/mlarson/L2/corsika/GeoCalibDetectorStatus_2012.56063_V1.i3.gz /gpfs/group/dfc13/default/dasha/mlarson/L2/corsika_all/Level2_IC86.2012_corsika.011808.$run_pad.i3.bz2  -o /gpfs/group/dfc13/default/dasha/mlarson/output_1/corsika/Pm_corsika_011808.$run_pad.part$subpart --sk $num_skip --ne $num_events 

end=`date +%s`
runtime=$((end-start))
echo 'time :'$runtime
