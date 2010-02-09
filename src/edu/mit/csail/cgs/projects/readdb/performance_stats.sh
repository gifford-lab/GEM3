#!/bin/bash

# sigmav7 chromlengths
chromstr="chr1:1-206082:+\nchr2:1-820529:+\nchr3:1-319572:+\nchr4:1-1489403:+\nchr5:1-589031:+\nchr6:1-253883:+\nchr7:1-1068014:+\nchr8:1-567694:+\nchr9:1-421853:+\nchr10:1-746106:+\nchr11:1-677346:+\nchr12:1-1095997:+\nchr13:1-928045:+\nchr14:1-774828:+\nchr15:1-1082448:+\nchr16:1-905116:+"

echo "retrieving all hits"
time echo -e "${chromstr}" | java edu.mit.csail.cgs.projects.readdb.Query --quiet --align 1102
time echo -e "${chromstr}" | java edu.mit.csail.cgs.projects.readdb.Query --quiet --align 1101
time echo -e "${chromstr}" | java edu.mit.csail.cgs.projects.readdb.Query --quiet --align 1082

echo "retrieving all hits and weights"
time echo -e "${chromstr}" | java edu.mit.csail.cgs.projects.readdb.Query --quiet --weights --align 1081
time echo -e "${chromstr}" | java edu.mit.csail.cgs.projects.readdb.Query --quiet --weights --align 1080
time echo -e "${chromstr}" | java edu.mit.csail.cgs.projects.readdb.Query --quiet --weights --align 1079

echo "retrieving histogram all chroms"
time echo -e "${chromstr}" | java edu.mit.csail.cgs.projects.readdb.Query --quiet --align 1077 --histogram 10
time echo -e "${chromstr}" | java edu.mit.csail.cgs.projects.readdb.Query --quiet --align 1076 --histogram 10
time echo -e "${chromstr}" | java edu.mit.csail.cgs.projects.readdb.Query --quiet --align 1075 --histogram 10

echo "random positions 100"
time ./random_positions.pl 100 | java edu.mit.csail.cgs.projects.readdb.Query --quiet --weights --align 1074
time ./random_positions.pl 100 | java edu.mit.csail.cgs.projects.readdb.Query --quiet --weights --align 1073
time ./random_positions.pl 100 | java edu.mit.csail.cgs.projects.readdb.Query --quiet --weights --align 1057

echo "random positions 10000"
time ./random_positions.pl 10000 | java edu.mit.csail.cgs.projects.readdb.Query --quiet --weights --align 1071
time ./random_positions.pl 10000 | java edu.mit.csail.cgs.projects.readdb.Query --quiet --weights --align 1070
time ./random_positions.pl 10000 | java edu.mit.csail.cgs.projects.readdb.Query --quiet --weights --align 1069

echo "simultaneous positions 10000 x 4"
for i in `seq 1052 1055` ; do
    time ./random_positions.pl 10000 | java edu.mit.csail.cgs.projects.readdb.Query --quiet --weights --align $i &
done

echo "simultaneous positions 10000 x 10"
for i in `seq 1052 1061` ; do
    time ./random_positions.pl 10000 | java edu.mit.csail.cgs.projects.readdb.Query --quiet --weights --align $i &
done


