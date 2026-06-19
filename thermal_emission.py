#!/usr/bin/env python3
"""
thermal_emission.py
-------------------
Sub-module for evaluating the direct thermal detectability of Planet Nine in the
3-4 micron band (WISE W1), replicating the low-temperature atmosphere models
of Fortney et al. (2016).

Assumes a methane-depleted stratosphere where methane has frozen/condensed out 
(abundance <= 10^-6), which opens a spectral window and increases the flux
significantly above blackbody expectations.
"""

import os
import argparse
import numpy as np
import matplotlib.pyplot as plt

# Constants
F0_W1 = 309.54 # Zero-point flux of WISE W1 (3.4 um) in Janskys
WISE_LIMIT_MAG = 16.5 # 5-sigma WISE All-Sky sensitivity limit in W1
WISE_LIMIT_FLUX = 10**6 * F0_W1 * 10**(-WISE_LIMIT_MAG / 2.5) # ~77.7 uJy

UNWISE_LIMIT_MAG = 17.0 # 5-sigma unWISE / NEOWISE coadd limit
UNWISE_LIMIT_FLUX = 10**6 * F0_W1 * 10**(-UNWISE_LIMIT_MAG / 2.5) # ~49.0 uJy


def compute_temperatures(d_hel, t_int=40.0, albedo=0.3):
    """
    Computes equilibrium temperature and effective temperature at distance d_hel (AU).
    T_eff = (T_int^4 + T_eq^4)^1/4
    T_eq = 278.3 * (1 - A)^0.25 / sqrt(d_hel)
    """
    t_eq = 278.3 * (1.0 - albedo)**0.25 / np.sqrt(d_hel)
    t_eff = (t_int**4 + t_eq**4)**0.25
    return t_eq, t_eff


def compute_w1_magnitude(d_hel, t_eff, methane_depleted=True):
    """
    Computes the apparent magnitude of Planet Nine in the WISE W1 (3.4 um) band
    using the Fortney et al. (2016) atmosphere models.
    
    Reference points from Fortney et al. 2016:
    For a methane-depleted atmosphere:
      - T_eff = 40 K, d = 700 AU -> W1 mag = 17.5
      - T_eff = 50 K, d = 700 AU -> W1 mag = 16.1
    We fit this scaling with:
      mag = mag_ref + 5 * log10(d / d_ref) - 15 * log10(T_eff / T_ref)
      
    For a methane-rich atmosphere (methane does not freeze out, e.g. Neptune-like):
      The flux is absorbed heavily, making the planet 6 magnitudes fainter (mag_ref ~ 23.5).
    """
    d_ref = 700.0
    t_ref = 40.0
    
    if methane_depleted:
        m_ref = 17.5
    else:
        m_ref = 23.5 # Methane absorbs 3-4 um flux heavily
        
    # Apparent magnitude scales:
    # 5 * log10(d / d_ref) due to geometric flux dilution (1/d^2 flux -> 5 log10(d) magnitude)
    # -15 * log10(T_eff / T_ref) matches Fortney's thermal flux temperature scaling (mag 17.5 to 16.1 from 40K to 50K)
    mag = m_ref + 5.0 * np.log10(d_hel / d_ref) - 15.0 * np.log10(t_eff / t_ref)
    return mag


def magnitude_to_flux_ujy(mag):
    """
    Converts apparent magnitude to flux density in microJanskys (uJy).
    """
    flux_jy = F0_W1 * 10**(-mag / 2.5)
    return flux_jy * 1e6


def run_thermal_model(a, e, t_int, albedo, output_dir):
    """
    Computes flux and magnitude as a function of true anomaly and plots curves.
    """
    os.makedirs(output_dir, exist_ok=True)
    
    true_anomaly = np.linspace(0, 360, 360) # Degrees
    true_anom_rad = np.radians(true_anomaly)
    
    # Orbit distance
    d_hel = a * (1.0 - e**2) / (1.0 + e * np.cos(true_anom_rad))
    
    # Calculate values for methane-depleted case
    t_eq, t_eff = compute_temperatures(d_hel, t_int, albedo)
    mag_depleted = compute_w1_magnitude(d_hel, t_eff, methane_depleted=True)
    flux_depleted = magnitude_to_flux_ujy(mag_depleted)
    
    # Calculate values for methane-rich case (standard atmosphere)
    mag_rich = compute_w1_magnitude(d_hel, t_eff, methane_depleted=False)
    flux_rich = magnitude_to_flux_ujy(mag_rich)
    
    # Diagnostic outputs
    idx_peri = 0 # Perihelion (f=0)
    idx_aph = 180 # Aphelion (f=180)
    
    print("\n--- Thermal Detection Diagnostic Summary ---")
    print(f"Heliocentric distance range: {d_hel[idx_peri]:.1f} AU (perihelion) to {d_hel[idx_aph]:.1f} AU (aphelion)")
    print(f"Effective temperature range: {t_eff[idx_peri]:.2f} K (perihelion) to {t_eff[idx_aph]:.2f} K (aphelion)")
    
    print("\n[Methane-Depleted Atmosphere Model (Fortney et al. 2016)]")
    print(f"  Perihelion: W1 mag = {mag_depleted[idx_peri]:.2f}, Flux = {flux_depleted[idx_peri]:.2f} uJy")
    print(f"  Aphelion:   W1 mag = {mag_depleted[idx_aph]:.2f}, Flux = {flux_depleted[idx_aph]:.2f} uJy")
    
    print("\n[Standard Methane-Rich Atmosphere Model]")
    print(f"  Perihelion: W1 mag = {mag_rich[idx_peri]:.2f}, Flux = {flux_rich[idx_peri]:.2f} uJy")
    print(f"  Aphelion:   W1 mag = {mag_rich[idx_aph]:.2f}, Flux = {flux_rich[idx_aph]:.2f} uJy")
    
    # Check detectability
    det_wise_peri = flux_depleted[idx_peri] >= WISE_LIMIT_FLUX
    det_wise_aph = flux_depleted[idx_aph] >= WISE_LIMIT_FLUX
    det_unwise_aph = flux_depleted[idx_aph] >= UNWISE_LIMIT_FLUX
    
    print("\n--- WISE / unWISE Detectability Vetting ---")
    print(f"  Can detect at perihelion (WISE)?      {'YES' if det_wise_peri else 'NO'}")
    print(f"  Can detect at aphelion (WISE)?        {'YES' if det_wise_aph else 'NO'}")
    print(f"  Can detect at aphelion (unWISE)?      {'YES' if det_unwise_aph else 'NO'}")
    
    # Plotting
    fig, axes = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
    
    # Subplot 1: Heliocentric Distance & Temperature
    ax1 = axes[0]
    color = 'tab:blue'
    ax1.set_ylabel('Heliocentric Distance (AU)', color=color)
    line1 = ax1.plot(true_anomaly, d_hel, color=color, linewidth=2, label='Heliocentric Distance')
    ax1.tick_params(axis='y', labelcolor=color)
    ax1.set_title("Planet Nine Orbital Phase & Thermal Emission Properties")
    
    ax2 = ax1.twinx()
    color = 'tab:red'
    ax2.set_ylabel('Effective Temperature (K)', color=color)
    line2 = ax2.plot(true_anomaly, t_eff, color=color, linestyle='--', linewidth=2, label='Effective Temp')
    ax2.tick_params(axis='y', labelcolor=color)
    
    lines = line1 + line2
    labels = [l.get_label() for l in lines]
    ax1.legend(lines, labels, loc='upper center')
    ax1.grid(True, linestyle=':', alpha=0.5)
    
    # Subplot 2: Flux density & Detection limits
    ax3 = axes[1]
    ax3.plot(true_anomaly, flux_depleted, 'r-', linewidth=2.5, label='Methane-Depleted Atmosphere')
    ax3.plot(true_anomaly, flux_rich, 'g--', linewidth=2.0, label='Methane-Rich (Neptune-like)')
    
    # Draw WISE / unWISE limits
    ax3.axhline(WISE_LIMIT_FLUX, color='black', linestyle=':', linewidth=1.5, label='WISE 5-sigma Limit (~77.7 uJy)')
    ax3.axhline(UNWISE_LIMIT_FLUX, color='grey', linestyle='-.', linewidth=1.5, label='unWISE 5-sigma Limit (~49.0 uJy)')
    
    ax3.set_xlabel('True Anomaly (degrees)')
    ax3.set_ylabel('Expected 3.4 um Flux Density (uJy)')
    ax3.set_yscale('log')
    ax3.set_ylim(1e-2, 1e3)
    ax3.set_xlim(0, 360)
    ax3.grid(True, which="both", linestyle=':', alpha=0.5)
    ax3.legend(loc='lower left')
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "thermal_emission_flux.png"), dpi=300)
    plt.close()
    
    # Save the second plot of Magnitude vs True Anomaly
    plt.figure(figsize=(8, 5))
    plt.plot(true_anomaly, mag_depleted, 'r-', linewidth=2.5, label='Methane-Depleted')
    plt.plot(true_anomaly, mag_rich, 'g--', linewidth=2.0, label='Methane-Rich')
    plt.axhline(WISE_LIMIT_MAG, color='black', linestyle=':', label='WISE limit (16.5)')
    plt.axhline(UNWISE_LIMIT_MAG, color='grey', linestyle='-.', label='unWISE limit (17.0)')
    plt.gca().invert_yaxis() # Magnitudes increase downwards (fainter)
    plt.xlabel('True Anomaly (degrees)')
    plt.ylabel('WISE W1 Apparent Magnitude')
    plt.title('Planet Nine WISE W1 Magnitude vs. Orbital Position')
    plt.xlim(0, 360)
    plt.grid(True, linestyle=':', alpha=0.5)
    plt.legend(loc='upper right')
    plt.savefig(os.path.join(output_dir, "thermal_emission_magnitude.png"), dpi=300)
    plt.close()
    
    print(f"Thermal emission figures saved to '{output_dir}'.")
    
    # Write diagnostic data to a text file
    with open(os.path.join(output_dir, "thermal_diagnostics.txt"), "w") as f:
        f.write("--- Thermal Detection Diagnostic Summary ---\n")
        f.write(f"a_9 = {a} AU, e_9 = {e}, T_int = {t_int} K, albedo = {albedo}\n")
        f.write(f"Heliocentric distance range: {d_hel[idx_peri]:.2f} AU (perihelion) to {d_hel[idx_aph]:.2f} AU (aphelion)\n")
        f.write(f"Effective temperature range: {t_eff[idx_peri]:.2f} K (perihelion) to {t_eff[idx_aph]:.2f} K (aphelion)\n\n")
        f.write("[Methane-Depleted Atmosphere Model]\n")
        f.write(f"  Perihelion: W1 mag = {mag_depleted[idx_peri]:.3f}, Flux = {flux_depleted[idx_peri]:.3f} uJy\n")
        f.write(f"  Aphelion:   W1 mag = {mag_depleted[idx_aph]:.3f}, Flux = {flux_depleted[idx_aph]:.3f} uJy\n\n")
        f.write("[Methane-Rich Atmosphere Model]\n")
        f.write(f"  Perihelion: W1 mag = {mag_rich[idx_peri]:.3f}, Flux = {flux_rich[idx_peri]:.3f} uJy\n")
        f.write(f"  Aphelion:   W1 mag = {mag_rich[idx_aph]:.3f}, Flux = {flux_rich[idx_aph]:.3f} uJy\n")


def main():
    parser = argparse.ArgumentParser(description="Planet Nine Methane-Depleted Thermal Emission Model")
    parser.add_argument("--a", type=float, default=550.0, help="Semi-major axis (AU)")
    parser.add_argument("--e", type=float, default=0.3, help="Eccentricity")
    parser.add_argument("--tint", type=float, default=40.0, help="Internal temperature (K)")
    parser.add_argument("--albedo", type=float, default=0.3, help="Geometric albedo")
    parser.add_argument("--output_dir", type=str, default="figures", help="Directory to save figures")
    args = parser.parse_args()
    
    run_model = run_thermal_model(args.a, args.e, args.tint, args.albedo, args.output_dir)


if __name__ == "__main__":
    main()
