#!/bin/bash 
#
# (0) Compile executables
make all
HOMEDIR=$(pwd)
#
# (1) Go to directory and link files
OUTPUT_DESTINATION=EXAMPLE/OUTPUT
mkdir -p ${OUTPUT_DESTINATION}
#
cd ${OUTPUT_DESTINATION}
ln -sf ${HOMEDIR}/bin/psphoton .
ln -sf ${HOMEDIR}/bin/post .
ln -sf ${HOMEDIR}/src/plot.py .
ln -sf ${HOMEDIR}/data/iasp91.smo .
cp ${HOMEDIR}/do.photon .
#
# (2) Run the phonon code with params from do file.
RUNTIME=15 # in seconds (run for several hours for real science)
gtimeout ${RUNTIME} bash do.photon
#
#
#
# (3) Postprocess the desired output file and write as ASCII tmp file.
TMPFILE=tmp.txt
./post << EOF
out.photon_zcore
${TMPFILE}
EOF
#
#
#
# (4) Read tmp file and plot array.
python plot.py << EOF
${TMPFILE}
EOF
#
#
#
# (5) Go back to original directory
cd ${HOMEDIR} 
