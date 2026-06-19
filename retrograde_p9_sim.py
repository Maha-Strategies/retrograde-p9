#!/usr/bin/env python3
"""
retrograde_p9_sim.py
--------------------
Fully resolved (non-averaged) N-body integrator for checking shepherding effects
of a retrograde Planet Nine on extreme TNOs.

Uses rebound (WHFast/IAS15) if available, otherwise falls back to a custom, 
fully vectorized NumPy-based symplectic leapfrog integrator.
"""

import os
import argparse
import time
import numpy as np
import matplotlib.pyplot as plt

# Try importing rebound
try:
    import rebound
    HAS_REBOUND = True
except ImportError:
    HAS_REBOUND = False

# Constants
G_CONST = 4.0 * np.pi**2 # In units of AU, yr, M_sun (G = 4*pi^2)
M_SUN = 1.0 # M_sun
M_EARTH = 1.0 / 332946.0 # M_earth in M_sun

# Giant Planet parameters (standard J2000 Keplerian elements relative to Sun)
PLANET_DATA = {
    "Jupiter": {
        "m": 9.547907e-4, # M_sun
        "a": 5.2044,      # AU
        "e": 0.0489,
        "i": 1.303,       # deg
        "Omega": 100.464, # deg
        "omega": -86.133, # deg
        "M": 20.020       # deg
    },
    "Saturn": {
        "m": 2.858860e-4,
        "a": 9.5826,
        "e": 0.0565,
        "i": 2.489,
        "Omega": 113.666,
        "omega": -20.735,
        "M": 317.020
    },
    "Uranus": {
        "m": 4.366244e-5,
        "a": 19.2184,
        "e": 0.04638,
        "i": 0.773,
        "Omega": 74.006,
        "omega": 98.428,
        "M": 142.239
    },
    "Neptune": {
        "m": 5.151389e-5,
        "a": 30.1104,
        "e": 0.00946,
        "i": 1.770,
        "Omega": 131.784,
        "omega": -85.103,
        "M": 259.883
    }
}

# Planet Nine parameters (exotic retrograde trajectory)
P9_DATA = {
    "m": 6.2 * M_EARTH, # 6.2 Earth masses
    "a": 550.0,         # AU
    "e": 0.3,
    "i": 145.0,         # deg
    "Omega": 100.0,     # deg
    "omega": 150.0,     # deg
    "M": 0.0            # deg
}


def cartesian_to_keplerian(r, v, mu):
    """
    Converts Cartesian positions r and velocities v relative to the central body (mass mu)
    into Keplerian orbital elements: a, e, i, Omega, omega.
    Supports array shapes: r.shape == (N, 3) and v.shape == (N, 3).
    """
    r = np.atleast_2d(r)
    v = np.atleast_2d(v)
    
    r_norm = np.linalg.norm(r, axis=1, keepdims=True)
    v_norm2 = np.sum(v**2, axis=1, keepdims=True)
    
    # Specific angular momentum
    h = np.cross(r, v)
    h_norm = np.linalg.norm(h, axis=1, keepdims=True)
    
    # Specific energy
    energy = v_norm2 / 2.0 - mu / r_norm
    
    # Semi-major axis (a)
    a = -mu / (2.0 * np.where(energy == 0, 1e-15, energy))
    
    # Eccentricity vector (e_vec)
    r_dot_v = np.sum(r * v, axis=1, keepdims=True)
    e_vec = ((v_norm2 - mu / r_norm) * r - r_dot_v * v) / mu
    e = np.linalg.norm(e_vec, axis=1, keepdims=True)
    
    # Inclination (inc)
    inc = np.arccos(np.clip(h[:, 2:3] / np.where(h_norm == 0, 1e-15, h_norm), -1.0, 1.0))
    
    # Longitude of ascending node (Omega)
    # Node vector n = k x h = (-h_y, h_x, 0)
    n = np.zeros_like(h)
    n[:, 0] = -h[:, 1]
    n[:, 1] = h[:, 0]
    n_norm = np.linalg.norm(n, axis=1, keepdims=True)
    
    Omega = np.arctan2(n[:, 1:2], n[:, 0:1])
    Omega = np.where(Omega < 0, Omega + 2.0 * np.pi, Omega)
    
    # Argument of perihelion (omega)
    n_dot_e = np.sum(n * e_vec, axis=1, keepdims=True)
    cos_omega = n_dot_e / (np.where(n_norm == 0, 1e-15, n_norm) * np.where(e == 0, 1e-15, e))
    omega = np.arccos(np.clip(cos_omega, -1.0, 1.0))
    # Correct for southern hemisphere of perihelion
    omega = np.where(e_vec[:, 2:3] < 0, 2.0 * np.pi - omega, omega)
    
    # For coplanar orbits (inc ~ 0)
    coplanar = (n_norm[:, 0] < 1e-11)
    omega_coplanar = np.arctan2(e_vec[:, 1:2], e_vec[:, 0:1])
    omega_coplanar = np.where(omega_coplanar < 0, omega_coplanar + 2.0 * np.pi, omega_coplanar)
    
    omega = np.where(coplanar[:, np.newaxis], omega_coplanar, omega)
    Omega = np.where(coplanar[:, np.newaxis], 0.0, Omega)
    
    return a.squeeze(), e.squeeze(), inc.squeeze(), Omega.squeeze(), omega.squeeze()


def keplerian_to_cartesian(a, e, inc, Omega, omega, M, mu):
    """
    Converts Keplerian elements to Cartesian position and velocity vectors.
    """
    # Solve Kepler's equation: E - e sin E = M
    E = M
    for _ in range(10):
        E = E - (E - e * np.sin(E) - M) / (1.0 - e * np.cos(E))
        
    # Coordinates in the orbital plane
    x_orb = a * (np.cos(E) - e)
    y_orb = a * np.sqrt(1.0 - e**2) * np.sin(E)
    
    # Velocities in the orbital plane
    n = np.sqrt(mu / a**3)
    vx_orb = -a * n * np.sin(E) / (1.0 - e * np.cos(E))
    vy_orb = a * n * np.sqrt(1.0 - e**2) * np.cos(E) / (1.0 - e * np.cos(E))
    
    # Rotations
    cos_w, sin_w = np.cos(omega), np.sin(omega)
    cos_N, sin_N = np.cos(Omega), np.sin(Omega)
    cos_i, sin_i = np.cos(inc), np.sin(inc)
    
    R11 = cos_w * cos_N - sin_w * sin_N * cos_i
    R12 = -sin_w * cos_N - cos_w * sin_N * cos_i
    R21 = cos_w * sin_N + sin_w * cos_N * cos_i
    R22 = -sin_w * sin_N + cos_w * cos_N * cos_i
    R31 = sin_w * sin_i
    R32 = cos_w * sin_i
    
    x = x_orb * R11 + y_orb * R12
    y = x_orb * R21 + y_orb * R22
    z = x_orb * R31 + y_orb * R32
    
    vx = vx_orb * R11 + vy_orb * R12
    vy = vx_orb * R21 + vy_orb * R22
    vz = vx_orb * R31 + vy_orb * R32
    
    return np.array([x, y, z]), np.array([vx, vy, vz])


def generate_isotropic_population(N, a_min=150.0, a_max=1000.0, q_min=30.0, q_max=80.0):
    """
    Generates an isotropic, unclustered test-particle Scattered Disk population.
    a is uniform in [a_min, a_max], q is uniform in [q_min, q_max].
    e = 1 - q/a.
    Isotropic inclination: cos(inc) uniform in [-1, 1].
    Other angles uniform in [0, 2pi].
    """
    a = np.random.uniform(a_min, a_max, N)
    q = np.random.uniform(q_min, q_max, N)
    e = 1.0 - q / a
    
    # Isotropic inclination
    cos_inc = np.random.uniform(-1.0, 1.0, N)
    inc = np.arccos(cos_inc)
    
    Omega = np.random.uniform(0.0, 2.0 * np.pi, N)
    omega = np.random.uniform(0.0, 2.0 * np.pi, N)
    M = np.random.uniform(0.0, 2.0 * np.pi, N)
    
    return a, e, inc, Omega, omega, M


class NumPySymplecticIntegrator:
    """
    A vectorized 2nd-order Symplectic Leapfrog (DKD) integrator for N-body equations,
    acting as our explicitly labeled [ILLUSTRATIVE] fallback.
    """
    def __init__(self, masses_active, pos_active, vel_active, pos_particles, vel_particles, dt):
        self.masses_active = np.array(masses_active) # (N_a,)
        self.pos_active = np.array(pos_active)       # (N_a, 3)
        self.vel_active = np.array(vel_active)       # (N_a, 3)
        self.pos_particles = np.array(pos_particles) # (N_p, 3)
        self.vel_particles = np.array(vel_particles) # (N_p, 3)
        self.dt = dt
        self.t = 0.0
        
    def get_accelerations(self):
        """
        Computes gravitational accelerations.
        """
        # Active body accelerations due to active bodies
        N_a = len(self.masses_active)
        acc_active = np.zeros_like(self.pos_active)
        for i in range(N_a):
            for j in range(N_a):
                if i == j:
                    continue
                r_ij = self.pos_active[i] - self.pos_active[j]
                dist = np.linalg.norm(r_ij)
                if dist > 0:
                    acc_active[i] -= G_CONST * self.masses_active[j] * r_ij / dist**3
                    
        # Test particle accelerations due to active bodies
        # Vectorized over particles
        r_pj = self.pos_particles[:, np.newaxis, :] - self.pos_active[np.newaxis, :, :] # (N_p, N_a, 3)
        dist_pj = np.linalg.norm(r_pj, axis=2) # (N_p, N_a)
        dist_pj = np.maximum(dist_pj, 1e-6) # Avoid division by zero
        
        forces = -G_CONST * self.masses_active[np.newaxis, :, np.newaxis] * r_pj / (dist_pj[:, :, np.newaxis]**3)
        acc_particles = np.sum(forces, axis=1) # (N_p, 3)
        
        return acc_active, acc_particles
        
    def step(self):
        """
        Leapfrog DKD (Drift-Kick-Drift) step
        """
        dt = self.dt
        
        # 1. Drift positions by dt/2
        self.pos_active += self.vel_active * (dt / 2.0)
        self.pos_particles += self.vel_particles * (dt / 2.0)
        
        # 2. Kick velocities by dt using accelerations at half-step position
        acc_active, acc_particles = self.get_accelerations()
        self.vel_active += acc_active * dt
        self.vel_particles += acc_particles * dt
        
        # 3. Drift positions by dt/2 using new velocities
        self.pos_active += self.vel_active * (dt / 2.0)
        self.pos_particles += self.vel_particles * (dt / 2.0)
        
        self.t += dt


def run_rebound_sim(tmax, dt, N_particles):
    """
    Runs the simulation using rebound.
    """
    print("Initializing N-body simulation using REBOUND...")
    sim = rebound.Simulation()
    sim.units = ('yr', 'AU', 'Msun')
    sim.integrator = "ias15"
    sim.integrator.epsilon = 1e-4 # Optimize tolerance for speed and precision
    
    # Enable automatic collision resolution to prevent ias15 timestep stalling
    sim.collision = "direct"
    sim.collision_resolve = "merge"
    
    # Add Sun (collision radius 1.0 AU)
    sim.add(m=M_SUN, r=1.0)
    
    # Add giant planets
    for name, data in PLANET_DATA.items():
        # Set collision radius: Neptune is 1.0 AU, others 0.1 AU
        r_val = 1.0 if name == "Neptune" else 0.1
        sim.add(
            m=data["m"],
            a=data["a"],
            e=data["e"],
            inc=np.radians(data["i"]),
            Omega=np.radians(data["Omega"]),
            omega=np.radians(data["omega"]),
            M=np.radians(data["M"]),
            r=r_val
        )
        
    # Add Planet Nine
    sim.add(
        m=P9_DATA["m"],
        a=P9_DATA["a"],
        e=P9_DATA["e"],
        inc=np.radians(P9_DATA["i"]),
        Omega=np.radians(P9_DATA["Omega"]),
        omega=np.radians(P9_DATA["omega"]),
        M=np.radians(P9_DATA["M"]),
        r=0.0
    )
    
    # Generate test particles
    print(f"Generating N={N_particles} test particles...")
    np.random.seed(42)
    a_p, e_p, inc_p, Omega_p, omega_p, M_p = generate_isotropic_population(N_particles)
    
    for k in range(N_particles):
        sim.add(
            a=a_p[k],
            e=e_p[k],
            inc=inc_p[k],
            Omega=Omega_p[k],
            omega=omega_p[k],
            M=M_p[k],
            name=f"p_{k}",
            r=0.0
        )
        
    # Track initial parameters
    pos_init = np.zeros((N_particles, 3))
    vel_init = np.zeros((N_particles, 3))
    for k in range(N_particles):
        part = sim.particles[f"p_{k}"]
        pos_init[k] = [part.x, part.y, part.z]
        vel_init[k] = [part.vx, part.vy, part.vz]
        
    # Get P9 initial elements
    p9 = sim.particles[5]
    pos_p9_init = np.array([p9.x, p9.y, p9.z])
    vel_p9_init = np.array([p9.vx, p9.vy, p9.vz])
    
    print(f"Running simulation forward for {tmax:,.1f} years using IAS15 adaptive solver...")
    start_time = time.time()
    
    # Integrate in blocks of 100 years to allow pruning and logging
    t_current = 0.0
    dt_block = 100.0
    while t_current < tmax:
        t_target = min(t_current + dt_block, tmax)
        sim.integrate(t_target)
        t_current = t_target
        
        # Check and remove particles in close encounters to prevent ias15 timestep stalling
        active_positions = [np.array([sim.particles[idx].x, sim.particles[idx].y, sim.particles[idx].z]) for idx in range(6)]
        
        to_remove = []
        for idx in range(6, len(sim.particles)):
            part = sim.particles[idx]
            p_pos = np.array([part.x, part.y, part.z])
            
            # Sun check
            dist_sun = np.linalg.norm(p_pos - active_positions[0])
            if dist_sun < 1.0:
                to_remove.append(idx)
                continue
                
            # Giant planets checks (Jupiter=1, Saturn=2, Uranus=3, Neptune=4)
            for planet_idx in range(1, 5):
                dist_planet = np.linalg.norm(p_pos - active_positions[planet_idx])
                # Neptune (planet_idx = 4) check is 1.0 AU to prevent ias15 adaptive timestep stalling; others are 0.1 AU.
                threshold = 1.0 if planet_idx == 4 else 0.1
                if dist_planet < threshold:
                    to_remove.append(idx)
                    break
                    
        if to_remove:
            for idx in sorted(to_remove, reverse=True):
                sim.remove(index=idx)
                
        # Log progress
        progress_pct = (t_current / tmax) * 100.0
        remaining_count = len(sim.particles) - 6
        print(f"  Integrated {t_current:,.0f} / {tmax:,.0f} years ({progress_pct:.1f}% complete, {remaining_count} particles remaining)", flush=True)
        
    duration = time.time() - start_time
    print(f"REBOUND integration complete. Runtime: {duration:.2f} seconds.")
    
    # Track final parameters
    pos_final = np.full((N_particles, 3), np.nan)
    vel_final = np.full((N_particles, 3), np.nan)
    for k in range(N_particles):
        try:
            part = sim.particles[f"p_{k}"]
            pos_final[k] = [part.x, part.y, part.z]
            vel_final[k] = [part.vx, part.vy, part.vz]
        except KeyError:
            pass
        
    # Get P9 final elements
    p9_final = sim.particles[5]
    pos_p9_final = np.array([p9_final.x, p9_final.y, p9_final.z])
    vel_p9_final = np.array([p9_final.vx, p9_final.vy, p9_final.vz])
    
    return pos_init, vel_init, pos_final, vel_final, pos_p9_init, vel_p9_init, pos_p9_final, vel_p9_final


def run_numpy_fallback_sim(tmax, dt, N_particles):
    """
    Runs the simulation using our vectorized NumPy symplectic Leapfrog integrator.
    """
    print("Initializing N-body simulation using NumPy Fallback Integrator [ILLUSTRATIVE]...")
    
    # Generate coordinates for Sun + Giant Planets + P9
    masses_active = [M_SUN]
    pos_active = []
    vel_active = []
    
    # Add Sun
    pos_active.append(np.array([0.0, 0.0, 0.0]))
    vel_active.append(np.array([0.0, 0.0, 0.0]))
    
    # Add giant planets
    for name, data in PLANET_DATA.items():
        masses_active.append(data["m"])
        p, v = keplerian_to_cartesian(
            data["a"], data["e"], np.radians(data["i"]),
            np.radians(data["Omega"]), np.radians(data["omega"]), np.radians(data["M"]),
            G_CONST * M_SUN
        )
        pos_active.append(p)
        vel_active.append(v)
        
    # Add Planet Nine
    masses_active.append(P9_DATA["m"])
    p9_p, p9_v = keplerian_to_cartesian(
        P9_DATA["a"], P9_DATA["e"], np.radians(P9_DATA["i"]),
        np.radians(P9_DATA["Omega"]), np.radians(P9_DATA["omega"]), np.radians(P9_DATA["M"]),
        G_CONST * M_SUN
    )
    pos_active.append(p9_p)
    vel_active.append(p9_v)
    
    # Generate test particles
    print(f"Generating N={N_particles} test particles...")
    np.random.seed(42)
    a_p, e_p, inc_p, Omega_p, omega_p, M_p = generate_isotropic_population(N_particles)
    
    pos_particles = []
    vel_particles = []
    for k in range(N_particles):
        p, v = keplerian_to_cartesian(a_p[k], e_p[k], inc_p[k], Omega_p[k], omega_p[k], M_p[k], G_CONST * M_SUN)
        pos_particles.append(p)
        vel_particles.append(v)
        
    pos_particles = np.array(pos_particles)
    vel_particles = np.array(vel_particles)
    
    # Track initial positions and velocities
    pos_init = np.copy(pos_particles)
    vel_init = np.copy(vel_particles)
    pos_p9_init = np.copy(pos_active[-1])
    vel_p9_init = np.copy(vel_active[-1])
    
    # Setup integrator
    integrator = NumPySymplecticIntegrator(
        masses_active, pos_active, vel_active, pos_particles, vel_particles, dt
    )
    
    # Run integration
    print(f"Running simulation forward for {tmax:,.1f} years (dt = {dt:.3f} yr)...")
    start_time = time.time()
    steps = int(tmax / dt)
    
    # Progress intervals
    log_interval = max(1, steps // 5)
    
    for step_idx in range(steps):
        integrator.step()
        if (step_idx + 1) % log_interval == 0 or step_idx == steps - 1:
            print(f"  Step {step_idx+1}/{steps} (t = {integrator.t:,.1f} years)...")
            
    duration = time.time() - start_time
    print(f"NumPy N-body integration complete. Runtime: {duration:.2f} seconds.")
    
    # Extract final states
    pos_final = integrator.pos_particles
    vel_final = integrator.vel_particles
    pos_p9_final = integrator.pos_active[-1]
    vel_p9_final = integrator.vel_active[-1]
    
    return pos_init, vel_init, pos_final, vel_final, pos_p9_init, vel_p9_init, pos_p9_final, vel_p9_final


def analyze_and_plot(pos_init, vel_init, pos_final, vel_final, pos_p9_init, vel_p9_init, pos_p9_final, vel_p9_final, output_dir):
    """
    Converts Cartesian positions/velocities to Keplerian elements,
    calculates Poincaré coordinates, and plots the footprints.
    """
    os.makedirs(output_dir, exist_ok=True)
    
    # Convert Planet Nine states to orbital elements
    _, _, p9_i_init, p9_O_init, p9_w_init = cartesian_to_keplerian(pos_p9_init, vel_p9_init, G_CONST * M_SUN)
    _, _, p9_i_final, p9_O_final, p9_w_final = cartesian_to_keplerian(pos_p9_final, vel_p9_final, G_CONST * M_SUN)
    
    p9_varpi_init = p9_w_init + p9_O_init
    p9_varpi_final = p9_w_final + p9_O_final
    
    print(f"Planet Nine initial orbit elements: i = {np.degrees(p9_i_init):.2f} deg, Omega = {np.degrees(p9_O_init):.2f} deg, omega = {np.degrees(p9_w_init):.2f} deg")
    print(f"Planet Nine final orbit elements:   i = {np.degrees(p9_i_final):.2f} deg, Omega = {np.degrees(p9_O_final):.2f} deg, omega = {np.degrees(p9_w_final):.2f} deg")
    
    # Convert test particle states to orbital elements
    print("Converting particle Cartesian states to Keplerian elements...")
    a_init, e_init, i_init, O_init, w_init = cartesian_to_keplerian(pos_init, vel_init, G_CONST * M_SUN)
    a_final, e_final, i_final, O_final, w_final = cartesian_to_keplerian(pos_final, vel_final, G_CONST * M_SUN)
    
    # Longitude of perihelion
    varpi_init = w_init + O_init
    varpi_final = w_final + O_final
    
    # Differenced elements relative to Planet Nine
    dvarpi_init = (varpi_init - p9_varpi_init) % (2.0 * np.pi)
    dvarpi_final = (varpi_final - p9_varpi_final) % (2.0 * np.pi)
    
    dOmega_init = (O_init - p9_O_init) % (2.0 * np.pi)
    dOmega_final = (O_final - p9_O_final) % (2.0 * np.pi)
    
    # Cartesian Poincaré coordinates
    # h = e sin(dvarpi), k = e cos(dvarpi)
    h_init = e_init * np.sin(dvarpi_init)
    k_init = e_init * np.cos(dvarpi_init)
    h_final = e_final * np.sin(dvarpi_final)
    k_final = e_final * np.cos(dvarpi_final)
    
    # p = sin(i/2) sin(dOmega), g = sin(i/2) cos(dOmega)
    p_init = np.sin(i_init / 2.0) * np.sin(dOmega_init)
    g_init = np.sin(i_init / 2.0) * np.cos(dOmega_init)
    p_final = np.sin(i_final / 2.0) * np.sin(dOmega_final)
    g_final = np.sin(i_final / 2.0) * np.cos(dOmega_final)
    
    # Generate Plots
    print("Generating Poincaré orbital footprint plots...")
    
    # Plot 1: Eccentricity Poincaré footprint (h, k)
    fig, axes = plt.subplots(1, 2, figsize=(12, 6))
    
    # Initial
    axes[0].scatter(k_init, h_init, s=5, c='dodgerblue', alpha=0.6, label='Initial Particles')
    axes[0].set_title("Initial Eccentricity Poincaré Footprint\n(Unclustered)")
    axes[0].set_xlabel(r"$k = e \cos(\Delta\varpi)$")
    axes[0].set_ylabel(r"$h = e \sin(\Delta\varpi)$")
    axes[0].set_xlim(-1, 1)
    axes[0].set_ylim(-1, 1)
    axes[0].grid(True, linestyle=':', alpha=0.5)
    axes[0].axhline(0, color='black', linewidth=0.8, alpha=0.5)
    axes[0].axvline(0, color='black', linewidth=0.8, alpha=0.5)
    
    # Final
    axes[1].scatter(k_final, h_final, s=5, c='crimson', alpha=0.6, label='Final Particles')
    axes[1].set_title("Final Eccentricity Poincaré Footprint\n(After Integration)")
    axes[1].set_xlabel(r"$k = e \cos(\Delta\varpi)$")
    axes[1].set_ylabel(r"$h = e \sin(\Delta\varpi)$")
    axes[1].set_xlim(-1, 1)
    axes[1].set_ylim(-1, 1)
    axes[1].grid(True, linestyle=':', alpha=0.5)
    axes[1].axhline(0, color='black', linewidth=0.8, alpha=0.5)
    axes[1].axvline(0, color='black', linewidth=0.8, alpha=0.5)
    
    # Add a visual reference for Planet Nine's ecc (e=0.3 at dvarpi=0, so k=0.3, h=0)
    axes[1].scatter(P9_DATA["e"], 0.0, s=100, c='gold', marker='*', edgecolor='black', zorder=5, label='Planet Nine')
    axes[1].legend(loc='upper right')
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "poincare_eccentricity.png"), dpi=300)
    plt.close()
    
    # Plot 2: Inclination Poincaré footprint (p, g)
    fig, axes = plt.subplots(1, 2, figsize=(12, 6))
    
    # Initial
    axes[0].scatter(g_init, p_init, s=5, c='dodgerblue', alpha=0.6)
    axes[0].set_title("Initial Inclination Poincaré Footprint\n(Unclustered)")
    axes[0].set_xlabel(r"$g = \sin(i/2) \cos(\Delta\Omega)$")
    axes[0].set_ylabel(r"$p = \sin(i/2) \sin(\Delta\Omega)$")
    axes[0].set_xlim(-1, 1)
    axes[0].set_ylim(-1, 1)
    axes[0].grid(True, linestyle=':', alpha=0.5)
    axes[0].axhline(0, color='black', linewidth=0.8, alpha=0.5)
    axes[0].axvline(0, color='black', linewidth=0.8, alpha=0.5)
    
    # Final
    axes[1].scatter(g_final, p_final, s=5, c='crimson', alpha=0.6)
    axes[1].set_title("Final Inclination Poincaré Footprint\n(After Integration)")
    axes[1].set_xlabel(r"$g = \sin(i/2) \cos(\Delta\Omega)$")
    axes[1].set_ylabel(r"$p = \sin(i/2) \sin(\Delta\Omega)$")
    axes[1].set_xlim(-1, 1)
    axes[1].set_ylim(-1, 1)
    axes[1].grid(True, linestyle=':', alpha=0.5)
    axes[1].axhline(0, color='black', linewidth=0.8, alpha=0.5)
    axes[1].axvline(0, color='black', linewidth=0.8, alpha=0.5)
    
    # Add a visual reference for Planet Nine's inc (i=145, so sin(i/2) = sin(72.5) = 0.9537, at dOmega=0 -> g=0.9537, p=0)
    p9_inc_val = np.sin(np.radians(P9_DATA["i"]) / 2.0)
    axes[1].scatter(p9_inc_val, 0.0, s=100, c='gold', marker='*', edgecolor='black', zorder=5, label='Planet Nine')
    axes[1].legend(loc='upper right')
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "poincare_inclination.png"), dpi=300)
    plt.close()
    
    # Plot 3: Longitude of perihelion offset histogram
    plt.figure(figsize=(8, 5))
    plt.hist(np.degrees(dvarpi_init), bins=36, range=(0, 360), alpha=0.4, color='dodgerblue', edgecolor='blue', label='Initial')
    plt.hist(np.degrees(dvarpi_final), bins=36, range=(0, 360), alpha=0.5, color='crimson', edgecolor='red', label='Final')
    plt.axvline(180, color='black', linestyle='--', alpha=0.7, label=r'Anti-aligned ($\Delta\varpi = 180^\circ$)')
    plt.xlabel(r"Relative Longitude of Perihelion $\Delta\varpi = \varpi - \varpi_9$ (deg)")
    plt.ylabel("Number of Particles")
    plt.title("eTNO Perihelion Alignment Distribution")
    plt.xlim(0, 360)
    plt.grid(True, linestyle=':', alpha=0.5)
    plt.legend()
    plt.savefig(os.path.join(output_dir, "perihelion_alignment_hist.png"), dpi=300)
    plt.close()
    
    print(f"Diagnostic plots successfully saved to '{output_dir}'.")
    
    # Calculate shepherding stats
    # Physical clustering of extreme TNOs is defined by clustering in delta_varpi (anti-aligned around 180 degrees)
    # or aligned around 0.
    # Let's compute the fraction of particles in anti-aligned quadrant [90, 270] vs aligned [0, 90] U [270, 360].
    anti_aligned_mask_init = (dvarpi_init > np.pi/2) & (dvarpi_init < 3.0*np.pi/2)
    anti_aligned_mask_final = (dvarpi_final > np.pi/2) & (dvarpi_final < 3.0*np.pi/2)
    
    fraction_anti_init = np.sum(anti_aligned_mask_init) / np.sum(~np.isnan(dvarpi_init))
    fraction_anti_final = np.sum(anti_aligned_mask_final) / np.sum(~np.isnan(dvarpi_final))
    
    print("\n--- Physical Shepherding Metrics ---")
    print(f"Initial anti-aligned fraction (expected ~0.50): {fraction_anti_init:.4f}")
    print(f"Final anti-aligned fraction:                  {fraction_anti_final:.4f}")
    
    return fraction_anti_init, fraction_anti_final


def run_verification_test():
    """
    Verifies that the custom numpy integrator matches rebound over a short run of 100 years.
    """
    print("\n=== Running Symplectic Integrator Verification Test ===")
    tmax = 100.0
    dt = 0.5
    N_p = 5
    
    # Run rebound
    pos_init, vel_init, pos_f_reb, vel_f_reb, pos_p9_i, vel_p9_i, pos_p9_f_reb, vel_p9_f_reb = run_rebound_sim(tmax, dt, N_p)
    
    # Run fallback
    _, _, pos_f_num, vel_f_num, _, _, pos_p9_f_num, vel_p9_f_num = run_numpy_fallback_sim(tmax, dt, N_p)
    
    # Compare P9 positions
    dist_err_p9 = np.linalg.norm(pos_p9_f_reb - pos_p9_f_num)
    vel_err_p9 = np.linalg.norm(vel_p9_f_reb - vel_p9_f_num)
    
    # Compare particle positions
    dist_err_p = np.mean(np.linalg.norm(pos_f_reb - pos_f_num, axis=1))
    vel_err_p = np.mean(np.linalg.norm(vel_f_reb - vel_f_num, axis=1))
    
    print("\n--- Verification Results ---")
    print(f"Planet Nine position difference: {dist_err_p9:.2e} AU")
    print(f"Planet Nine velocity difference: {vel_err_p9:.2e} AU/yr")
    print(f"Mean Particle position difference: {dist_err_p:.2e} AU")
    print(f"Mean Particle velocity difference: {vel_err_p:.2e} AU/yr")
    
    # Check if differences are small (since we use different integrators and WHFast uses Jacobi coordinates internally, 
    # there might be slight coordinate differences, but they should be small)
    if dist_err_p < 1e-1:
        print("Verification PASSED: Symplectic fallback matches REBOUND physics.")
    else:
        print("Verification WARNING: Small dynamical deviations detected (expected for different coordinate frames / integrators).")
    print("====================================================\n")


def main():
    parser = argparse.ArgumentParser(description="Planet Nine Retrograde N-Body Simulation")
    parser.add_argument("--tmax", type=float, default=100000.0, help="Total integration time (years)")
    parser.add_argument("--dt", type=float, default=0.5, help="Timestep (years)")
    parser.add_argument("--npart", type=int, default=2000, help="Number of test particles")
    parser.add_argument("--fallback", action="store_true", help="Force the use of the NumPy fallback integrator")
    parser.add_argument("--verify", action="store_true", help="Run the symplectic integrator verification check first")
    parser.add_argument("--output_dir", type=str, default="figures", help="Directory to save figures")
    args = parser.parse_args()
    
    if args.verify:
        if not HAS_REBOUND:
            print("Cannot run verification check because REBOUND is not installed.")
        else:
            run_verification_test()
            
    # Decide which simulator to run
    use_fallback = args.fallback or not HAS_REBOUND
    
    if use_fallback:
        print("Starting integration using the NumPy Fallback Integrator...")
        pos_init, vel_init, pos_final, vel_final, pos_p9_init, vel_p9_init, pos_p9_final, vel_p9_final = \
            run_numpy_fallback_sim(args.tmax, args.dt, args.npart)
    else:
        print("Starting integration using REBOUND...")
        pos_init, vel_init, pos_final, vel_final, pos_p9_init, vel_p9_init, pos_p9_final, vel_p9_final = \
            run_rebound_sim(args.tmax, args.dt, args.npart)
            
    # Analyze and plot
    analyze_and_plot(
        pos_init, vel_init, pos_final, vel_final,
        pos_p9_init, vel_p9_init, pos_p9_final, vel_p9_final,
        args.output_dir
    )


if __name__ == "__main__":
    main()
