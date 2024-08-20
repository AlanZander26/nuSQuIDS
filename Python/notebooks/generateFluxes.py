import matplotlib.pyplot as plt
import timeit
import sys

module_path = '../src/'

sys.path.append(module_path)
from implement_TeVSGT import *

start = timeit.default_timer()

model_type = "ADD"
params_ranges = [[0.2, 0.2, 1], [0, 0, 1]]
neutrino_type = "antineutrino" 
L = 12742 # km
E_min = 1 # GeV
E_max = 1e5 # GeV
E_range = [E_min, E_max]
medium = "earth"
initial_flux_ratios = [0, 1, 0]
NormalOrdering = True

n_mesh = 110
data = implement_model(model_type, params_ranges, E_range, [medium, [-1, 0, n_mesh]], initial_flux_ratios, neutrino_type, NormalOrdering, Nen=n_mesh)

energy = data[0, 0, 0, :, 0] # in GeV
cos = np.linspace(-1, 0, n_mesh)

flux_ratio_grid = data[:, 0, 0, :, 2].T

cosX, EY = np.meshgrid(cos, energy)

fig, ax = plt.subplots(1, 1, figsize=(12,6))


ax.set_yscale('log')
ax.set_xlabel("cos theta_z")
ax.set_ylabel("E (in GeV)")
ax.set_ylim([E_min, E_max])


plt.pcolormesh(cosX, EY, flux_ratio_grid, cmap = "RdYlBu_r", vmin=0, vmax=1)

fig.savefig('TeVSGT_Earth_angle.pdf',bbox_inches='tight', dpi=150)

stop = timeit.default_timer()

print('Time: ', stop - start)