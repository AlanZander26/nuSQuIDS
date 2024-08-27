import numpy as np
import os
import subprocess
import h5py

def implement_model(model_type, params_ranges, E_range, medium_list, initial_flux_ratios, neutrino_type="neutrino", NormalOrdering=True, Nen=200, Nen_grid=100):
    
    str_sed = '1c\\#define USE_' + model_type.upper()  
    nuSQuIDS_PATH = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    TEVSGT_PATH = nuSQuIDS_PATH + "/TeVSGT"
    
    subprocess.run(
        ["sed", "-i", str_sed, "main.cpp"],
        cwd=TEVSGT_PATH,
        text=True
    )
    
    # Convert parameters for the models into their correct types
    for sublist in params_ranges:
        sublist[0] = float(sublist[0])  # a or N
        sublist[1] = float(sublist[1])  # m0
        if len(sublist) == 3:
            sublist[2] = float(sublist[2])  # mu (if applicable)
        
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
    for sublist in params_ranges:
        for param in sublist:
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
    
    
    # Print the output data
    hdf5_file = os.path.expanduser(f"{TEVSGT_PATH}/outputfile.hdf5")
    with h5py.File(hdf5_file, 'r') as hdf_file:
        # Function to recursively print the structure of the HDF5 file
        def print_structure(name, obj):
            obj_type = "Group" if isinstance(obj, h5py.Group) else "Dataset"
            print(f"{name}: {obj_type}")
            
            # If it's a dataset, print its shape and dtype
            if isinstance(obj, h5py.Dataset):
                print(f"  Shape: {obj.shape}, DataType: {obj.dtype}")

        # Iterate over all objects in the file and print their structure
        hdf_file.visititems(print_structure)
    
    return