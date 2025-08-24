#TOV SOLVER
#Author: SUDIPTA HENSH
#EMAIL: SUDIPTAHENSH2009@GMAIL.COM

import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import PchipInterpolator
import matplotlib.pyplot as plt
import os

# Constants in CGS units
G = 6.67430e-8       # Gravitational constant
c = 2.99792458e10    # Speed of light
Msol = 1.98847e33    # Solar mass
kmtocm = 1e5         # km to cm conversion

def tovrhs(r, y):
    p, m = y
    if p <= 1e-10:
        return [0, 0]
    
    rho = eos(p)
    
    if r <= 1e-8:
        return [0, 0]
        
    if rho <= 0.0:
        print("density is negative", r)
        return [0, 0]
        
    if m <= 0.0:
        print("Mass is negative", r)
        return [0, 0]
    
    metric_term = 1 - 2*G*m/(c**2*r)
    # Stricter check near horizon formation
    if metric_term <= 1e-14:
        print(f"Terminating: Metric became non-physical")
        return [0, 0]
    
    # TOV equations
    dpdr = -(G*(rho + p/c**2)*(m + 4*np.pi*r**3*p/c**2))/(r**2 * metric_term)
    dmdr = 4*np.pi*r**2 * rho
    
    return [dpdr, dmdr]

def read_eos_table(filename):
    data = np.loadtxt(filename)
    pressure = data[:, 3]
    density = data[:, 2]
    
    # Convert to log-log space
    log_p = np.log10(pressure + 1e-10)
    log_rho = np.log10(density + 1e-10)
    
    # Use PCHIP interpolation (monotonic)
    interp = PchipInterpolator(log_p, log_rho)
    
    return lambda p: 10**interp(np.log10(p + 1e-10)) if p > 0 else 0

def solve_tov(p_center, max_radius=250):
    r0 = 1e-8
    rho_center = eos(p_center)
    m0 = (4/3)*np.pi*r0**3 * rho_center
    
    # Adaptive surface threshold
    surface_threshold = max(1e5, 1e-10 * p_center)  # Higher threshold for high pressures
    
    def surface_event(r, y):
        return y[0] - surface_threshold
    
    # Stricter solver parameters for high pressures
    if p_center > 1e35:
        rtol_val = 1e-8 # relative tolerance
        atol_val = 1e-8 # absolute tolerance
        max_step_val = 1e4  # Smaller steps (1 km)
        method = 'Radau'
    else:
        rtol_val = 1e-8
        atol_val = 1e-8
        max_step_val = 1e5  # Larger steps (10 km)
        method = 'RK45'
        
    sol = solve_ivp(tovrhs, [r0, max_radius*kmtocm], [p_center, m0],
                   events=surface_event,
                   rtol=rtol_val, atol=atol_val,
                   method=method,
                   max_step=max_step_val)
    
    if sol.t_events[0].size > 0:
        return (sol.t_events[0][0]/kmtocm, sol.y_events[0][0][1]/Msol)
    return (np.nan, np.nan)

# Main execution
eos_file = './eos_tables/QHC21D.lorene'
eos_table = np.loadtxt(eos_file)
pressure_data = eos_table[:, 3]
density_data = eos_table[:, 2]

min_pressure = 1.2e33 #minimum pressure in dyn/cm^2
max_pressure = 1e37 #maximum pressure in dyn/cm^2
max_radius=30 # above this integration won't continue and a solution won't be found.
#usually this limit is good enough for M-R relation of a TOV star

print(f"Pressure range: {min_pressure:.3e} - {max_pressure:.3e} dyn/cm²")

eos = read_eos_table(eos_file)

# Create filename based on EoS file name
base_name = os.path.splitext(os.path.basename(eos_file))[0]
output_dir = './out_data'
# Ensure directory exists
os.makedirs(output_dir, exist_ok=True)
output_filename = os.path.join(output_dir, base_name + '_MR_relation.txt')

# Create pressure array with more points in high-pressure region
pc_low = np.logspace(np.log10(min_pressure), np.log10(1e35), 180) # less points when pressure is below 1e35 dyn/cm^2
pc_high = np.logspace(np.log10(1e35), np.log10(max_pressure), 180)  # more points when pressure is above 1e35 dyn/cm^2
pc = np.unique(np.concatenate([pc_low, pc_high]))

radii, masses = [], []
successful_pressures = []

for p in pc:
    r, m = solve_tov(p, max_radius=30)
    if not np.isnan(r):
        radii.append(r)
        masses.append(m)
        successful_pressures.append(p)
        print(f"p_center = {p:.3e} -> R = {r:.3f} km, M = {m:.3f} Msun")
    else:
        print(f"Integration failed at p_center = {p:.3e}")

# Find maximum mass
if masses:
    max_mass_idx = np.argmax(masses)
    max_mass = masses[max_mass_idx]
    max_mass_radius = radii[max_mass_idx]
    max_mass_pressure = successful_pressures[max_mass_idx]
    print(f"\nMaximum mass: {max_mass:.3f} M☉ at R = {max_mass_radius:.3f} km")
    print(f"Central pressure for maximum mass: {max_mass_pressure:.3e} dyn/cm²")
    
    
# Save results to file
if radii:  # Only save if we have valid results
   
    # Prepare data - columns: central pressure, radius, mass
    data_to_save = np.column_stack((successful_pressures, radii, masses))
    
    # Create header with column descriptions and parameters
    header = (
        "Mass-Radius Relation for Neutron Stars\n"
        f"EoS file: {eos_file}\n"
        "Columns:\n"
        "1. Central Pressure (dyn/cm²)\n"
        "2. Radius (km)\n"
        "3. Mass (Solar Masses)\n"
        f"Min pressure: {min_pressure:.3e}, Max pressure: {max_pressure:.3e}\n"
        f"Number of points: {len(radii)}"
    )
    
    # Save to file
    np.savetxt(output_filename, data_to_save,
               fmt=['%.6e', '%.6f', '%.6f'],
               header=header,
               delimiter='\t')
    
    print(f"\nSaved mass-radius data to: {output_filename}")
    

# Plot
plt.figure(figsize=(10, 6))
plt.plot(radii, masses, 'o-', markersize=4)

# Highlight maximum mass
if masses:
    plt.plot(max_mass_radius, max_mass, 'ro', markersize=8, label=f'Max Mass: {max_mass:.3f} M☉')

# Plot settings
plt.xlabel('Radius (km)', fontsize=14)
plt.ylabel('Mass (M$_\odot$)', fontsize=14)
plt.title('Neutron Star Mass-Radius Relation', fontsize=16)
plt.grid(True, alpha=0.3)
plt.legend()
plt.tight_layout()

plt.show()
