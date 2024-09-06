import numpy as np
import os
import subprocess

def implement_model(model_type, params, E_range, medium_list, initial_flux_ratios, neutrino_type="neutrino", NormalOrdering=True, Nen=200, Nen_grid=100):
    
    str_sed = '1c\\#define USE_' + model_type.upper()  
    nuSQuIDS_PATH = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    TEVSGT_PATH = nuSQuIDS_PATH + "/TeVSGT"
    
    subprocess.run(
        ["sed", "-i", str_sed, "main.cpp"],
        cwd=TEVSGT_PATH,
        text=True
    )
    
    params = [float(p) for p in params]
        
    # Process energy range
    E_range = np.array(E_range, dtype=float)
    E_min, E_max = E_range
    
    # Get medium and medium parameter details
    medium = medium_list[0]
    medium_param_range = medium_list[-1]
    medium_param_min = float(medium_param_range[0])
    medium_param_max = float(medium_param_range[1])
    N_medium_param = int(medium_param_range[2])
    
    # Convert initial flux ratios
    initial_flux_ratios = np.array(initial_flux_ratios, dtype=float)
    
    # Convert NormalOrdering to string representation
    NormalOrdering = str(NormalOrdering).lower()
    
    # Construct the input string for the bash script
    input_str = ""
    for param in params:
        input_str += f"{param}\n"
    
    input_str += f"{E_min}\n{E_max}\n{medium}\n{medium_param_min}\n{medium_param_max}\n{N_medium_param}\n{neutrino_type}\n{NormalOrdering}\n{Nen_grid}\n{Nen}\n"
    
    for ratio in initial_flux_ratios:
        input_str += f"{ratio}\n"
    
    # Run the bash script and provide the input
    subprocess.run(
        ["bash", f"{TEVSGT_PATH}/run_commands.sh"],
        input=input_str,
        text=True
    )
    
    return