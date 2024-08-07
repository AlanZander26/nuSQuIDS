import numpy as np
import os
import subprocess


def implement_model(model_type, params_ranges, E_range, medium_list, initial_flux_ratios, neutrino_type="neutrino", NormalOrdering=True, Nen=1000, Nen_grid=200):
    
    str_sed = '1c\\#define USE_' + model_type.upper()  
    nuSQuIDS_PATH = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    TEVSGT_PATH = nuSQuIDS_PATH + "/TeVSGT"
    



    subprocess.run(
        ["sed", "-i", str_sed, "main.cpp"],
        cwd=TEVSGT_PATH,
        text=True
    )
    

    
    for sublist in params_ranges:
        sublist[0] = float(sublist[0])
        sublist[1] = float(sublist[1])
        sublist[2] = int(sublist[2])
        
    tuple_shape = ()
    for sublist in params_ranges[::-1]:
        tuple_shape += tuple([sublist[2]])
    tuple_shape += tuple([Nen, 3+1])
    
    
    params_ranges = [item for sublist in params_ranges for item in sublist]
    E_range = np.array(E_range, dtype=float)
    E_min, E_max = E_range
    medium = medium_list[0]
    medium_param_range = medium_list[-1]
    medium_param_min = float(medium_param_range[0])
    medium_param_max = float(medium_param_range[1])
    N_medium_param = int(medium_param_range[2])
    initial_flux_ratios = np.array(initial_flux_ratios, dtype=float)
    NormalOrdering = str(NormalOrdering).lower()
    
    # Construct the input string for the bash script
    input_str = ""
    for param_input in params_ranges:
        input_str += f"{param_input}\n"
    input_str += f"{E_min}\n{E_max}\n{medium}\n{medium_param_min}\n{medium_param_max}\n{N_medium_param}\n{neutrino_type}\n{NormalOrdering}\n{Nen_grid}\n{Nen}\n"
    for ratio in initial_flux_ratios:
        input_str += f"{ratio}\n"
    
    # Run the bash script and provide the input
    subprocess.run(
        ["bash", f"{TEVSGT_PATH}/run_commands.sh"],

    input=input_str,
    text=True
    )
    
    tuple_shape = tuple([N_medium_param]) + tuple_shape
    data_file = os.path.expanduser(f"{TEVSGT_PATH}/fluxes_flavor.txt")
    data = np.loadtxt(data_file).reshape(tuple_shape)
    return data