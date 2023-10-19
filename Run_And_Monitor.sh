#!/bin/bash

# Function to monitor and log max memory usage
monitor_memory() {
  local job_id=$1
  local log_file=$2
  local max_memory=0
  
  while :; do
    # Get memory usage for this job
    memory=$(sacct -j $job_id --format=MaxRSS --noheader | awk '{print $1}' | sed '/^$/d' | sed 's/K//')
  
    # Update max_memory if needed
    if [ "$memory" -gt "$max_memory" ]; then
      max_memory=$memory
    fi
    
    # Log the max memory usage so far
    echo "Max memory used so far: ${max_memory}K" >> $log_file
  
    sleep 5  # Sleep for 5 seconds (adjust as needed)
  done
}

if [ "$#" -lt 2 ]; then
  echo "Usage: $0 <command> <log_file_path>"
  exit 1
fi

# Your command to run
command_to_run=$1

# Your log file path
log_file_path=$2

# Create a temporary SLURM script
echo "#!/bin/bash" > temp_slurm_script.sh
echo "$command_to_run" >> temp_slurm_script.sh

# Submit job and capture the output
output=$(sbatch --wrap "bash temp_slurm_script.sh")

# Parse the output to get the SLURM job ID
job_id=$(echo $output | awk '{print $NF}')

# Log or use the job ID
echo "SLURM Job ID: $job_id"


# Wait for the job to complete
while [ "$(scontrol show job $job_id | grep 'JobState=' | awk -F= '{print $2}' | awk '{print $1}')" != "COMPLETED" ]; do
    sleep 5  # Wait for 5 seconds before checking again
done

# Start memory monitoring after the job is complete
monitor_memory $job_id $log_file_path

sleep 5
#Get node name to monitor process
node_name=$(squeue -j $job_id -o %N --noheader)


# Wait for the job to complete
srun -w "$node_name" -J "wait_$job_id" sleep 1


# Wait for the job to complete
#srun -w "$(squeue -j $job_id -o %N)" -J "wait_$job_id" sleep 1

# Terminate the memory monitoring function
kill $!

# Clean up
#rm temp_slurm_script.sh
