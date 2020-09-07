#!/bin/bash     

NUMBERS=$(seq 1 1)
for NUM in ${NUMBERS}
do
echo "Submitting: ${NUM}"

while IFS=$'\t' read P1
do 
JOB=`qsub << EOF

# -A Rarray.out
#PBS -q lisandro
# -m e
# -M lbenedetti@biologia.unipi.it
#PBS -j oe        
#PBS -N marine_heat_waves  
#PBS -l nodes=1:ppn=72 

module load R/3.5.1

Rscript ~/MHWs/code/HeatWaves.R $1 ${P1}
EOF
`
echo "JobID = ${JOB} for parameters ${P1} 
#submitted on `date`"
echo PBS_NODEFILE is $PBS_NODEFILE
done < ~/MHWs/data/mhw_parms_${NUM}.txt

done
echo "Script completed"            
