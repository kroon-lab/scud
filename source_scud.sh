#-------------------------------
# Add custom path to python path
#-------------------------------
#export myscud=/home/bate/scud
export PYTHONPATH=${PYTHONPATH}:${PWD}

#---------------------------
# source cctbx (from phenix)
#---------------------------

# Local phenix distribution should contain this file:
source /eiwit/progs/phenix/phenix/phenix_env.sh
export PHENIX_TRUST_OTHER_ENV=1

#-----------
# add to bin
#-----------

export PATH=${PATH}:${PWD}/bin

