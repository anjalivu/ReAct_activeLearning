from UC_Function import UCModel
from inputParams import *  # Your parameter file
import csv
import os

L_mean = L
w_f_mean = w_f
t_t_mean = t_t
t_FRP_mean = t_FRP
prestress_mean = prestress

N = 10 # Number of simulations to run
stdev_factor = .25

for job_id in range(N):

    # add randomness to parameters
    L = np.random.normal(L_mean, L_mean * stdev_factor) 
    w_f = np.random.normal(w_f_mean, w_f_mean * stdev_factor)  
    t_t = np.random.normal(t_t_mean, t_t_mean * stdev_factor)  
    t_FRP = np.random.normal(t_FRP_mean, t_FRP_mean * stdev_factor) 
    prestress = np.random.normal(prestress_mean, prestress_mean * stdev_factor)  

    # Create unique job name for each simulation
    job_name = f'Job-{job_id}'

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
    
    # Clean up all files except .odb
    cleanup_extensions = ['.dat', '.msg', '.sta', '.com', '.prt', '.sim', '.log', '.lck']
    
    for ext in cleanup_extensions:
        filename = job_name + ext
        if os.path.exists(filename):
            try:
                os.remove(filename)
                print(f"Removed {filename}")
            except Exception as e:
                print(f"Could not remove {filename}: {e}")
    
    print(f"Cleanup complete for {job_name}. ODB file retained.\n")