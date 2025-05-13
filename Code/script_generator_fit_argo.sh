#!/bin/bash
# SLURM template
read -r -d '' SLURM_TEMPLATE << EOM
#!/bin/bash -l
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=14
#SBATCH --partition=batch
#SBATCH -J {GridID}_{modelType}_{presLevel}
#SBATCH -o logs/{GridID}_{modelType}_{presLevel}.out
#SBATCH -e logs/{GridID}_{modelType}_{presLevel}.err
#SBATCH --time=240:00:00
#SBATCH --mem=6G

module load udunits/2.2.28 
module load proj/9.1.1/gnu-12.2.0
module load gdal/3.6.2/gnu-12.2.0 
module load geos/3.12.1 
module load java/19.0.1
module load llvm/16.0.1
module load R/4.3.0/gnu-12.2.0
export R_LIBS_USER=/home/saduakd/local/R4.3.0_libs.gnu:$R_LIBS_USER

module load gcc
export OMP_NUM_THREADS=14
export OMP_PLACES=cores
export OMP_PROC_BIND=close

srun -c 14 Rscript ../fit_argo_nig.R {GridID} {modelType} {presLevel}
EOM

GridIDs="374"
modelTypes="gauss_indep gauss_cor nig_indep nig_cor"
presLevels="1000"

for GridID in $GridIDs; do
    for presLevel in $presLevels; do
        for modelType in $modelTypes; do
            JOB_FILE="slurm_${GridID}_${modelType}_${presLevel}.sh"
            echo "${SLURM_TEMPLATE//\{GridID\}/$GridID}" | \
                sed "s/{modelType}/$modelType/g"| \
                sed "s/{presLevel}/$presLevel/g" > $JOB_FILE

            # Submit the job file to SLURM
            sbatch $JOB_FILE
        done
    done
done
