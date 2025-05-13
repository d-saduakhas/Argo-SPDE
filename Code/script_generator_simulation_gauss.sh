#!/bin/bash

# Define the ranges for your parameters
rhos=(2 1 1 -2)
rhos_eps=(-0.8 -0.4 -0.1 0 0.1 0.4 0.8)
n_obs=1000
n_rep=10

# Create logs directory if it doesn't exist
mkdir -p logs

# SLURM template
read -r -d '' SLURM_TEMPLATE << EOM
#!/bin/bash -l
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=28
#SBATCH --partition=batch
#SBATCH -J g_rho_{rho}_rhoeps_{rho_eps}_$SLURM_JOB_ID
#SBATCH -o logs/g_rho_{rho}_rhoeps_{rho_eps}_$SLURM_JOB_ID.out
#SBATCH -e logs/g_rho_{rho}_rhoeps_{rho_eps}_$SLURM_JOB_ID.err
#SBATCH --time=150:00:00
#SBATCH --mem=10G

module load udunits/2.2.28 
module load proj/9.1.1/gnu-12.2.0
module load gdal/3.6.2/gnu-12.2.0 
module load geos/3.12.1 
module load java/19.0.1
module load llvm/16.0.1
module load R/4.3.0/gnu-12.2.0
export R_LIBS_USER=/home/saduakd/local/R4.3.0_libs.gnu:$R_LIBS_USER
module load gcc
export OMP_NUM_THREADS=28
export OMP_PLACES=cores
export OMP_PROC_BIND=close

srun Rscript simulation_study_gauss.R {rho} {rho_eps} $n_obs $n_rep
EOM

# Iterate over the parameter combinations and create individual job scripts
for rho in "${rhos[@]}"; do
    for rho_eps in "${rhos_eps[@]}"; do
        JOB_FILE="slurm_rho_${rho}_rhoeps_${rho_eps}.sh"
        echo "${SLURM_TEMPLATE//\{rho\}/$rho}" | \
            sed "s/{rho_eps}/$rho_eps/g" > $JOB_FILE
        
        # Submit the job file to SLURM using sbatch
        sbatch $JOB_FILE
    done
done

echo "All simulations submitted."
