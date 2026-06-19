#!/usr/bin/env python3
"""
run_all.py
----------
Coordination script that runs the retrograde Planet Nine N-body simulation
and the methane-depleted thermal emission model in sequence, saving all
outputs to the figures directory.
"""

import os
import sys
import subprocess
import time

def main():
    print("==========================================================")
    print("   Planet Nine Retrograde Orbit N-Body Simulation Runner  ")
    print("==========================================================\n")
    
    start_time = time.time()
    
    # Create output directories
    output_dir = "figures"
    os.makedirs(output_dir, exist_ok=True)
    
    # We choose a default of 20,000 years to run the full proof-of-concept simulation
    # while showing the initial secular footprint evolution for 2,000 particles.
    tmax = 20000.0
    dt = 0.5
    npart = 2000
    
    print(f"--- Step 1: Running N-body integration (T_max = {tmax:,.0f} yr, N = {npart}) ---")
    sim_cmd = [
        sys.executable, "retrograde_p9_sim.py",
        "--tmax", str(tmax),
        "--dt", str(dt),
        "--npart", str(npart),
        "--output_dir", output_dir
    ]
    
    sim_start = time.time()
    try:
        subprocess.run(sim_cmd, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error executing N-body simulation: {e}")
        sys.exit(1)
    sim_duration = time.time() - sim_start
    print(f"N-body simulation step complete in {sim_duration:.2f} seconds.\n")
    
    # 2. Run the Thermal Emission Module
    print("--- Step 2: Running Thermal Emission Flux Calculations ---")
    thermal_cmd = [
        sys.executable, "thermal_emission.py",
        "--a", "550.0",
        "--e", "0.3",
        "--tint", "40.0",
        "--albedo", "0.3",
        "--output_dir", output_dir
    ]
    
    thermal_start = time.time()
    try:
        subprocess.run(thermal_cmd, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error executing thermal emission calculation: {e}")
        sys.exit(1)
    thermal_duration = time.time() - thermal_start
    print(f"Thermal emission calculation complete in {thermal_duration:.2f} seconds.\n")
    
    total_duration = time.time() - start_time
    print("==========================================================")
    print("   Execution Completed Successfully!   ")
    print(f"   Total Time: {total_duration:.2f} seconds")
    print(f"   All plots saved to the '{output_dir}' directory.")
    print("==========================================================\n")

if __name__ == "__main__":
    main()
