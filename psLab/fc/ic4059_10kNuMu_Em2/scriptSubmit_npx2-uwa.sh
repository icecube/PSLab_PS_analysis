#!/bin/sh

# 2008/09/09  CBF
# A modification of Submit_npx2-uwa for scripts,
# so that filename numbering matches script output file numberings

#   ./scriptSubmit_npx2-uwa  255  thisShellScript000255.sh

# A simple sh script for submitting one job to the npx2-uwa cluster.
# Simply log into npx2-uwa.icecube.wisc.edu, run any command as you normally 
#  would but prepend ./submit-npx2-uwa.sh and the command will be run on the
#  npx2-uwa cluster.
# This script will create directories to store your execution script, log files,
#  errors, and std output, so you need write permission.

# This script creates a script to be executed and another script to submit it.
# The execution script must be available *at time of job execution!*, which may
#  not be until much later. 
# The submission script can be discarded immediately.

  # Creating output directories
  mkdir -p uwa-execs uwa-logs uwa-out uwa-error

  # Creating execution script, do not delete until job has started!
  echo "#!/bin/bash" > uwa-execs/uwa-$1.sh
  echo "date" >> uwa-execs/uwa-$1.sh
  echo "hostname" >> uwa-execs/uwa-$1.sh
  echo "cd `pwd`" >> uwa-execs/uwa-$1.sh
  echo "$2" >> uwa-execs/uwa-$1.sh
  echo "date" >> uwa-execs/uwa-$1.sh

  chmod +x uwa-execs/uwa-$1.sh

  # Creating condor submission script (ClassAd)
  echo "Universe  = vanilla" > 2sub.sub
  echo "Executable = uwa-execs/uwa-$1.sh" >> 2sub.sub
  echo "Log = uwa-logs/uwa-$1.log" >> 2sub.sub
  echo "Output = uwa-out/uwa-$1.out" >> 2sub.sub
  echo "Error = uwa-error/uwa-$1.error" >> 2sub.sub
  echo "Notification = NEVER" >> 2sub.sub 
  echo "Queue" >> 2sub.sub
  condor_submit 2sub.sub


