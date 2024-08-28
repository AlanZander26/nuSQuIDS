import timeit
import sys

module_path = '../src/'

sys.path.append(module_path)
from implement_TeVSGT import *

start = timeit.default_timer()

# Changeable Parameters
model_type = sys.argv[1]
params_ranges = eval(sys.argv[2])
n_mesh = int(sys.argv[3])
n_final_mesh = int(sys.argv[4])
E_min = float(sys.argv[5])
E_max = float(sys.argv[6])
NormalOrdering = True

# Fixed Parameters
medium = "earth"
L = 12742 # km
initial_flux_ratios = [0, 1, 0]
neutrino_type = "antineutrino"
E_range = [E_min, E_max]


implement_model(model_type, params_ranges, E_range, [medium, [-1, 0, n_mesh]], initial_flux_ratios, neutrino_type, NormalOrdering, Nen=n_final_mesh, Nen_grid=n_mesh)

stop = timeit.default_timer()

print('Time: ', stop - start)