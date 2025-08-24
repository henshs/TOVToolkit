# TOVToolkit
Python-based solver for the Tolman–Oppenheimer–Volkoff equations with tabulated equations of state.

This repository contains a Python script to solve the Tolman-Oppenheimer-Volkoff (TOV) equations for neutron stars using a tabulated equation of state (EoS). It computes the mass-radius (M-R) relation for a range of central pressures and plots the results.

## Features

* Solves TOV equations for a given tabulated EoS in CGS units.
* Uses adaptive numerical integration with SciPy's `solve_ivp`.
* Supports high-precision integration for high central pressures.
* Reads EoS tables (compatible with LORENE-like formats, easy to adapt if format is not the same).
* Saves M-R relation to an output text file.
* Plots the mass-radius curve and highlights the maximum mass.

## Prerequisites

Make sure you have the following installed:

* Python 3.8+
* NumPy
* SciPy
* Matplotlib

Install dependencies with:

```bash
pip install numpy scipy matplotlib
```

## Input Data Format

The EoS file should be a plain text table with at least four columns. The script expects:

* Column 2: Mass density (g/cm³)
* Column 3: Pressure (dyn/cm²)
* Comments in the file should start with `#` (ignored).
* Adjust your table or the script as per your need.
* There are a bunch of EoS tables in LORENE format.

Example: `./eos_tables/eos_sly4.lorene`

## How to Run

1. Place your EoS file in the `./eos_tables` directory.
2. Update the `eos_file` path in the script to point to your EoS file.
3. Run the script:

   ```bash
   python tovtoolkit.py
   ```
4. The script will:

   * Print the mass-radius data to the console.
   * Save the results to `./out_data/<EoS_name>_MR_relation.txt`.
   * Plot the mass-radius relation.

## Output

* **Text file**: Contains central pressure, radius, and mass for each integration.
* **Plot**: Displays the neutron star mass-radius relation.

## Notes

* Adjust `min_pressure`, `max_pressure`, and `max_radius` in the script to explore different mass-radius ranges.
* The script uses monotonic PCHIP interpolation to avoid spurious oscillations in the density-pressure relation.
* Metric breakdown or negative densities trigger early termination.

## References

The script integrates the standard TOV equations:

Tolman, R.C. (1939). "Static Solutions of Einstein's Field Equations for Spheres of Fluid". *Physical Review*, 55(4), 364–373. [DOI:10.1103/PhysRev.55.364](https://doi.org/10.1103/PhysRev.55.364)

Oppenheimer, J.R., Volkoff, G.M. (1939). "On Massive Neutron Cores". *Physical Review*, 55(4), 374–381. [DOI:10.1103/PhysRev.55.374](https://doi.org/10.1103/PhysRev.55.374)

## License

This project is open-source. You are free to use, modify, and distribute it with attribution.

