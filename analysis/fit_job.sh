#! /bin/bash
#SBATCH -A dp004
#SBATCH -p cosma7
#SBATCH --job-name=fit_df
#SBATCH --output=out/pd.master.snap40.%a.%N.%j.o
#SBATCH --error=err/pd.master.snap40.%a.%N.%j.e
# #SBATCH --mail-type=ALL
# #SBATCH --mail-user=c.lovell@herts.ac.uk
#SBATCH -t 0-03:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
# #SBATCH --mem-per-cpu=3800
#SBATCH --array=0-5

# module purge
# module load gnu_comp/7.3.0
# module load hdf5/1.10.3


module load pythonconda3/4.5.4
source activate eagle

python fit_gsmf.py $SLURM_ARRAY_TASK_ID
# python fit_sfrf.py $SLURM_ARRAY_TASK_ID

echo "Job done, info follows..."
sacct -j $SLURM_JOBID --format=JobID,JobName,Partition,Elapsed,ExitCode,MaxRSS,CPUTime,SystemCPU,ReqMem
exit
