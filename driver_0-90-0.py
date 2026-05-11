from UC_Conditional import UC_conditional
import csv
import os
import numpy as np
from datetime import datetime
import time
import glob

# model geometry
L = 100.0

# Material properties
#FRP
E2_FRP = 7010.0
nu12_FRP = 0.314
G12_FRP = 4661.0
G13_FRP = G12_FRP
G23_FRP = 2540.0
rho_FRP = 1.2e-9
# TPU
rho_m = 1.0e-9
C10_m = 30.972
D1_m = 0.0

# tape
rho_t = 1.2e-9
E_t = 1.0
nu_t = 0.3

# FRP layup
t_t = 0.13

# out of plane pulling displacement, selected to achieve same rotation as a simulation that is known to work
s_alpha = np.sqrt(2)*8.0/12.7

# job setup
cpus = 1






N = 50
theta_range = [-90, 90]
phi_range = [-90, 90]
t_frp_range = [0.025, 0.060]
w_f_range = [6.0, 15.0]
prestress_range = [4.0, 8.0]
E1_FRP_range = [50000.0, 250000.0]

n_dim = 6

seed = 67
np.random.seed(seed)

bins = np.arange(N)
samples = np.zeros((N, n_dim))
for i in range(n_dim):
    np.random.shuffle(bins)
    # (bin + random[0,1)) / N ensures values stay in [0, 1)
    samples[:, i] = (bins) / N
noise = np.random.uniform(0, 1/N, size=samples.shape)
samples = samples + noise

for sample in samples:
    theta = 0.0
    phi = 90.0
    t_FRP = np.interp(sample[2], [0, 1], t_frp_range)
    w_f = np.interp(sample[3], [0, 1], w_f_range)
    prestress = np.interp(sample[4], [0, 1], prestress_range)
    E1_FRP = np.interp(sample[5], [0, 1], E1_FRP_range)

    uz_pull = w_f/np.sqrt(2)*s_alpha
    layup = [theta, phi, theta]
    meshSize = w_f/5.0





    job_name = f'data_gen_0-90-0'

    # Call the function with unique job name
    radius, has_inflection, pull_corners_needed = UC_conditional(L, w_f, E1_FRP, E2_FRP, nu12_FRP, G12_FRP, 
                                    G13_FRP, G23_FRP, rho_FRP, rho_m, C10_m, D1_m, 
                                    rho_t, E_t, nu_t, t_t, t_FRP, layup, meshSize, 
                                    prestress, uz_pull, cpus, job_name)

    # Check results
    print(f"Cylinder radius: {radius}")
    print(f"Is tristable: {not has_inflection}")

    # Define CSV file path
    csv_file = 'training_data_0-90-0.csv'

    # Check if file exists to determine if we need to write headers
    file_exists = os.path.isfile(csv_file)

    # Open CSV file in append mode and write results
    with open(csv_file, 'a', newline='') as f:
        writer = csv.writer(f)
        
        # Write header if file doesn't exist
        if not file_exists:
            writer.writerow(['theta', 'phi', 't_FRP', 'w_f', 'prestress', 'E1_FRP', 'radius', 'pull_corners_needed', 'has_inflection'])
        
        # Write data row
        writer.writerow([theta, phi, t_FRP, w_f, prestress, E1_FRP, radius, pull_corners_needed, has_inflection])

    print(f"Results saved to {csv_file}")
    
    # Cleanup simulation files
    for pattern in [f'{job_name}.*', f'{job_name}_pullCorners.*']:
        for filepath in glob.glob(pattern):
            if not filepath.endswith('.cae'):
                try:
                    os.remove(filepath)
                except:
                    pass