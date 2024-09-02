import h5py
import os
import re
import numpy as np

directory = os.path.dirname(os.path.abspath(__file__))

output_directory = os.path.join(directory, 'precalculatedFluxes')
os.makedirs(output_directory, exist_ok=True)

# Initialize lists to store zenith angles (CosTheta values)
zenith_angles = []

# Find a sample filename to derive the combined file name
sample_filename = None
for filename in os.listdir(directory):
    if filename.endswith('.hdf5') and 'CosTheta' in filename:
        sample_filename = filename
        break

# Derive the combined file name by removing the 'CosTheta_(angle)' part
if sample_filename:
    combined_file_name = re.sub(r'_CosTheta_-?\d+\.\d+', '', sample_filename)
    combined_file_path = os.path.join(output_directory, combined_file_name)
else:
    raise ValueError("No HDF5 files with 'CosTheta' found in the directory.")

# Create a new HDF5 file to store the combined data
with h5py.File(combined_file_path, 'w') as combined_file:
    
    # Variable to check if the energy range has been copied
    energy_range_copied = False

    # List all HDF5 files in the directory
    for filename in os.listdir(directory):
        if filename.endswith('.hdf5') and 'CosTheta' in filename and filename != combined_file_name:
            # Extract CosTheta value from the filename
            match = re.search(r'CosTheta_(-?\d+\.\d+)', filename)
            if match:
                costheta = float(match.group(1))
                group_name = f'costh_{costheta}'
                zenith_angles.append(costheta)
                
                # Open the current HDF5 file
                file_path = os.path.join(directory, filename)
                with h5py.File(file_path, 'r') as current_file:
                    
                    # Copy the energy_range dataset if not already done
                    if not energy_range_copied:
                        combined_file.create_dataset('energy_range', data=current_file['energies'][...])
                        energy_range_copied = True
                    
                    # Create a group in the combined file with the extracted CosTheta value
                    group = combined_file.create_group(group_name)
                    
                    # Copy all datasets and groups from the current file to the new group
                    def copy_attrs(obj, group):
                        for key, value in obj.attrs.items():
                            group.attrs[key] = value

                    def copy_datasets(name, node):
                        if isinstance(node, h5py.Dataset):
                            combined_file[group_name + '/' + name] = node[...]
                        elif isinstance(node, h5py.Group):
                            grp = combined_file.create_group(group_name + '/' + name)
                            copy_attrs(node, grp)
                    
                    current_file.visititems(copy_datasets)
                
                # After processing the file, delete it
                os.remove(file_path)
    
    # Now create the zenith_angles dataset with the collected CosTheta values
    combined_file.create_dataset('zenith_angles', data=np.array(zenith_angles))

print(f'All HDF5 files have been combined into {combined_file_path}, and the original files have been deleted.')
