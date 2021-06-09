#-------------------------------------------------------------
#script for running a single-CPU serial job via SLURM
#-------------------------------------------------------------
#
#SBATCH --job-name=sandpiper
#SBATCH --output=sandpiper_output
#
#Define the number of hours the job should run. 
#Maximum runtime is limited to 10 days, ie. 240 hours
#SBATCH --time=36:00:00
#
#Define the amount of RAM used by your job in GigaBytes
#SBATCH --mem=2G
#
#Send emails when a job starts, it is finished or it exits
#SBATCH --mail-user=raimundojulian.saonaurmeneta@ist.ac.at
#SBATCH --mail-type=ALL
#
#Pick whether you prefer requeue or not. 
# Options: --requeue or --no-requeue
#SBATCH --no-requeue
#
#Do not export the local environment to the compute nodes
#SBATCH --export=NONE
unset SLURM_EXPORT_ENV
#
#Make sure that they use a single thread
export OMP_NUM_THREADS=1
#
#load the respective software module you intend to use
module load rust
#
#
#run the respective binary through SLURM's srun
srun --cpu_bind=verbose  cluster.exe