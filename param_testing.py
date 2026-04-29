from UC_Function import UCModel
from inputParams import *  # Your parameter file
import csv
import os
import numpy as np
from scipy.stats import truncnorm
from datetime import datetime
import time

def get_truncated_normal(mean, sd, low, high):
    """
    Generate a truncated normal distribution.
    
    Parameters:
    mean: mean of the normal distribution
    sd: standard deviation of the normal distribution
    low: lower bound
    high: upper bound
    """
    a = (low - mean) / sd
    b = (high - mean) / sd
    return truncnorm(a, b, loc=mean, scale=sd)

# Store original mean values
L_mean = L
w_f_mean = w_f
t_t_mean = t_t
t_FRP_mean = t_FRP
prestress_mean = prestress

N = 20  # Number of simulations to run
stdev_factor = 0.25

for i in range(N):
    # Generate timestamp-based job ID (format: MMDD_HHMM_iteration)
    timestamp = datetime.now().strftime("%m%d_%H%M")
    job_id = f"{timestamp}_{i:03d}"  # e.g., "1225_1430_000"
    
    # Add randomness to parameters using truncated normal
    # This prevents values from going negative or too far from the mean
    
    # L parameter - truncated between 50% and 150% of mean
    L_dist = get_truncated_normal(
        mean=L_mean, 
        sd=L_mean * stdev_factor, 
        low=L_mean * 0.5, 
        high=L_mean * 1.5
    )
    L = L_dist.rvs()
    
    # w_f parameter - truncated between 50% and 150% of mean
    w_f_dist = get_truncated_normal(
        mean=w_f_mean, 
        sd=w_f_mean * stdev_factor, 
        low=w_f_mean * 0.5, 
        high=w_f_mean * 1.5
    )
    w_f = w_f_dist.rvs()
    
    # t_t parameter - truncated between 50% and 150% of mean
    t_t_dist = get_truncated_normal(
        mean=t_t_mean, 
        sd=t_t_mean * stdev_factor, 
        low=t_t_mean * 0.5, 
        high=t_t_mean * 1.5
    )
    t_t = t_t_dist.rvs()
    
    # t_FRP parameter - truncated between 50% and 150% of mean
    t_FRP_dist = get_truncated_normal(
        mean=t_FRP_mean, 
        sd=t_FRP_mean * stdev_factor, 
        low=t_FRP_mean * 0.5, 
        high=t_FRP_mean * 1.5
    )
    t_FRP = t_FRP_dist.rvs()
    
    # prestress parameter - truncated between 50% and 150% of mean
    prestress_dist = get_truncated_normal(
        mean=prestress_mean, 
        sd=prestress_mean * stdev_factor, 
        low=prestress_mean * 0.5, 
        high=prestress_mean * 1.5
    )
    prestress = prestress_dist.rvs()

    # Create unique job name for each simulation
    job_name = f'Job-{job_id}'
    
    print(f"\n{'='*60}")
    print(f"Starting simulation {i+1}/{N}")
    print(f"Job name: {job_name}")
    print(f"{'='*60}\n")

    # Call the function with unique job name
    radius, has_inflection = UCModel(L, w_f, E1_FRP, E2_FRP, nu12_FRP, G12_FRP, 
                                    G13_FRP, G23_FRP, rho_FRP, rho_m, C10_m, D1_m, 
                                    rho_t, E_t, nu_t, t_t, t_FRP, layup, meshSize, 
                                    prestress, uz_pull, cpus, job_name)

    # Check results
    print(f"Cylinder radius: {radius}")
    print(f"Is tristable: {not has_inflection}")
    
    # Export results to CSV
    csv_filename = "simulation_results.csv"
    file_exists = os.path.isfile(csv_filename)
    
    with open(csv_filename, 'a', newline='') as csvfile:
        fieldnames = ['Job_ID','L', 'w_f', 'E1_FRP', 'E2_FRP', 'nu12_FRP', 'G12_FRP', 
                    'G13_FRP', 'G23_FRP', 'rho_FRP', 'rho_m', 'C10_m', 'D1_m', 
                    'rho_t', 'E_t', 'nu_t', 't_t', 't_FRP', 'layup', 'meshSize', 
                    'prestress', 'uz_pull', 'cpus', 'Radius', 'Has_Inflection', 'Is_Tristable']
        
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        
        # Write header only if file is new
        if not file_exists:
            writer.writeheader()
        
        # Write results
        writer.writerow({
            'Job_ID': job_id,
            'L': L,
            'w_f': w_f,
            'E1_FRP': E1_FRP,
            'E2_FRP': E2_FRP,
            'nu12_FRP': nu12_FRP,
            'G12_FRP': G12_FRP,
            'G13_FRP': G13_FRP,
            'G23_FRP': G23_FRP,
            'rho_FRP': rho_FRP,
            'rho_m': rho_m,
            'C10_m': C10_m,
            'D1_m': D1_m,
            'rho_t': rho_t,
            'E_t': E_t,
            'nu_t': nu_t,
            't_t': t_t,
            't_FRP': t_FRP,
            'layup': layup,
            'meshSize': meshSize,
            'prestress': prestress,
            'uz_pull': uz_pull,
            'cpus': cpus,
            'Radius': radius,
            'Has_Inflection': has_inflection,
            'Is_Tristable': not has_inflection
        })
    
    print(f"Results saved to {csv_filename}")
    
    # Wait a moment to ensure all files are written
    time.sleep(2)
    
    # Clean up files immediately after this job completes
    cleanup_extensions = ['.dat', '.msg', '.sta', '.com', '.prt', '.sim', '.log', '.lck', '.023', '.pac', '.res', '.stt', '.mdl', '.abq']
    
    # Try to clean up with the job_name
    files_removed = 0
    for ext in cleanup_extensions:
        filename = job_name + ext
        if os.path.exists(filename):
            try:
                os.remove(filename)
                print(f"Removed {filename}")
                files_removed += 1
            except Exception as e:
                print(f"Could not remove {filename}: {e}")
    
    # Also check for files without the "Job-" prefix (in case Abaqus strips it)
    for ext in cleanup_extensions:
        filename = job_id + ext
        if os.path.exists(filename):
            try:
                os.remove(filename)
                print(f"Removed {filename}")
                files_removed += 1
            except Exception as e:
                print(f"Could not remove {filename}: {e}")
    
    # Check for any numbered job files (job-0, job-1, etc.) and remove them too
    for ext in cleanup_extensions:
        for j in range(10):  # Check job-0 through job-9
            filename = f'job-{j}{ext}'
            if os.path.exists(filename):
                try:
                    os.remove(filename)
                    print(f"Removed {filename}")
                    files_removed += 1
                except Exception as e:
                    print(f"Could not remove {filename}: {e}")
    
    print(f"Cleanup complete for {job_name}. Removed {files_removed} files. ODB file retained.")
    print(f"Completed job {i+1}/{N}\n")

print(f"\n{'='*60}")
print(f"All {N} simulations completed!")
print(f"{'='*60}\n")