#!/bin/sh
#PBS -W group_list=dfc13_collab
#PBS -l nodes=1:ppn=1
#PBS -l mem=10gb
#PBS -l walltime=4:00:00
#PBS -j oe
#PBS -o /storage/home/dup193/work/estes/output/test1.log
#PBS -d /storage/home/dup193/work/estes/new_estes/STV
#PBS -N nue
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
ENV_SHELL="/gpfs/group/dfc13/default/dasha/icerec_V05-01-05_TH/build/env-shell.sh"
SCRIPT="/storage/home/dup193/work/estes/new_estes/STV/Extract.py"
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



outdir=/storage/home/dup193/work/estes/new_estes/STV/

if [ ! -d "$DIRECTORY" ]; then
    mkdir -p ${outdir}
fi
cd ${outdir}
p1= $1
# p2= $2
echo $p1
# echo $p2

i3_run py2-v2 "$ENV_SHELL" 'python' "$SCRIPT" $1 #$2 

end=`date +%s`
runtime=$((end-start))
echo 'time :'$runtime
