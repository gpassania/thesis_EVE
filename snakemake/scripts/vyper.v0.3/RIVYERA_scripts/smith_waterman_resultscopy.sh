#!/bin/bash

############################################################################
# Generic Rivyera shell script to copy the smith-waterman results to
# a target server/folder, and create a "signal" file "*.done" when the 
# copying has completed.
# This is done from the Rivyera, because of empty-file problems when
# trying to copy the results data from the linux cluster node. 
############################################################################
# Michael Forster                              m.forster@ikmb.uni-kiel.de 
# v0.3.3   2 September 2014
############################################################################

rsafile=/home/mforster/.ssh/nukey

if [ $# -ne 3 ]
then
echo "You must give three parameters: <smithwaterman results file> <donefile> <target>"
echo " e.g. smith_waterman_resultscopy.sh /home/mforster/data/A503_smith-waterman_detailed.txt /home/mforster/data/A503_out.fpga.fa.done user@server.uni-kiel.de:/targetfolder"
exit
fi

############################################################################
###  smith-waterman detailed results file, "done" file, and
###  target server/folder location from command line arguments
############################################################################

detailedfile=${1}
donefile=${2}
targetlocation=${3}

############################################################################
###  copy the results accross to target server/folder
############################################################################

# wait 5 seconds in case of operating system sychronisation problems on Rivyera:
sleep 5

# copy smith-waterman detailed results file over to target server/folder:
##this waits for a password because it does not find the rsa key:
##rsync -av ${detailedfile} ${targetlocation}
rsync -av -e "ssh  -i $rsafile" ${detailedfile} ${targetlocation}

############################################################################
### create a flag-file and copy this accross to target server/folder
############################################################################

touch ${donefile}

# wait 15 seconds in case of operating system sychronisation problems on target server:
sleep 15

# copy flag-file over to target server/folder:
##rsync -av ${donefile} ${targetlocation}
rsync -av -e "ssh  -i $rsafile" ${donefile} ${targetlocation}

# wait 15 seconds in case of operating system sychronisation problems on target server:
sleep 15
