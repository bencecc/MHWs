#!/bin/bash     
while IFS=$'\t' read P1 P2

do 
JOB=`qsub << EOF

# -A Rarray.out
#PBS -q lisandro
# -m e
# -M lbenedetti@biologia.unipi.it
#PBS -j oe        
#PBS -N R_array_test  
#PBS -l nodes=1:ppn=72 

module load R/3.5.1

Rscript ~/MHWs/code/rand_phase_net.R $1 $2 ${P1}
EOF
`
echo "JobID = ${JOB} for parameters ${P1}
#submitted on `date`"
#echo PBS_NODEFILE is $PBS_NODEFILE
done < ~/MHWs/data/null_net_1.txt              
