from UC_Function import UCModel
from inputParams import *  # Your parameter file
import csv
import os

prestresses = [6.0, 6.01, 6.02]
for job_id in range(len(prestresses)):
    prestress = prestresses[job_id]
    # Call the function
    radius, has_inflection = UCModel(L, w_f, E1_FRP, E2_FRP, nu12_FRP, G12_FRP, 
                                    G13_FRP, G23_FRP, rho_FRP, rho_m, C10_m, D1_m, 
                                    rho_t, E_t, nu_t, t_t, t_FRP, layup, meshSize, 
                                    prestress, uz_pull, cpus, job_id)

    # Check results
    print(f"Cylinder radius: {radius}")
    print(f"Is tristable: {not has_inflection}")
    
    # Export results to CSV
    csv_filename = "simulation_results.csv"
    file_exists = os.path.isfile(csv_filename)
    
    with open(csv_filename, 'a', newline='') as csvfile:
        fieldnames = ['Job ID','L', 'w_f', 'E1_FRP', 'E2_FRP', 'nu12_FRP', 'G12_FRP', 
                    'G13_FRP', 'G23_FRP', 'rho_FRP', 'rho_m', 'C10_m', 'D1_m', 
                    'rho_t', 'E_t', 'nu_t', 't_t', 't_FRP', 'layup', 'meshSize', 
                    'prestress', 'uz_pull', 'cpus', 'Radius', 'Has_Inflection', 'Is_Tristable']
        
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        
        # Write header only if file is new
        if not file_exists:
            writer.writeheader()
        
        # Write results
        writer.writerow({
            'Job ID': job_id,
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
