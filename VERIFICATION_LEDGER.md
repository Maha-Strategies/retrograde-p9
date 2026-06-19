# Verification Ledger: Mathematical & Physical Formulations

This ledger records the mathematical steps, physical equations, and coordinate transformations used in the N-body integrations and thermal emission calculations for the retrograde Planet Nine configuration.

---

## 1. N-Body Equations of Motion

### Active Point-Mass Formulation (Our Model)
In our simulation (`retrograde_p9_sim.py`), we model the Sun, Jupiter, Saturn, Uranus, Neptune, and Planet Nine as fully resolved active point masses. The equations of motion for each active body $i$ are:

$$\frac{d^2 \mathbf{r}_i}{dt^2} = -G \sum_{j \neq i} m_j \frac{\mathbf{r}_i - \mathbf{r}_j}{|\mathbf{r}_i - \mathbf{r}_j|^3}$$

For each massless test particle $k$, the equation of motion is:

$$\frac{d^2 \mathbf{w}_k}{dt^2} = -G \sum_{j \in \text{active}} m_j \frac{\mathbf{w}_k - \mathbf{r}_j}{|\mathbf{w}_k - \mathbf{r}_j|^3}$$

Where:
* $G = 4\pi^2\text{ AU}^3\text{ yr}^{-2}\text{ M}_\odot^{-1}$ is the gravitational constant.
* $m_j$ is the mass of active body $j$.
* $\mathbf{r}_j$ is the barycentric position vector of active body $j$.
* $\mathbf{w}_k$ is the position vector of test particle $k$.

---

## 2. Comparison with Historical Secular Approximations

In historical literature (e.g., Batygin et al. 2019), many long-term dynamical analyses utilize the **phase-averaged secular approximation** (or quadrupolar wire approximation). Here we flag the critical structural differences between these two methodologies:

| Physical Phenomena | Active Point-Mass Calculations (Our Model) | Semi-Averaged Secular Equations (Batygin et al. 2019) |
| :--- | :--- | :--- |
| **Planetary Representation** | Active moving point masses with varying orbital phases. | Orbits of Jupiter, Saturn, Uranus, and Neptune are averaged (spread out) into static concentric rings ("wires" of mass). |
| **Gravitational Potential** | Fully time-dependent, resolved gravitational forces: $V(\mathbf{r}) = -G \sum \frac{m_j}{|\mathbf{r} - \mathbf{r}_j|}$. | Axisymmetric quadrupole-smoothed potential ($J_2$ smoothing of the central Sun-planet system). |
| **Resonant Dynamics** | Captures **Mean-Motion Resonances (MMRs)** (e.g. 3:1, 2:1 resonances) which can trap and shepherd eTNOs. | Filters out mean anomalies, completely removing all MMRs and resonant trapping mechanics. |
| **Close Encounters** | Accurately models **gravitational scattering** (close flybys of Neptune or other giant planets), which alters $a$ and $q$. | Particles can pass directly "through" the planetary wires without experiencing localized scattering, artificially conserving $a$ and $q$. |
| **Integration Timescales** | Uses REBOUND's adaptive, non-symplectic `ias15` solver. IAS15 manages step sizes dynamically to ensure numerical precision ($10^{-15}$ local error tolerance) during retrograde crossings, bypassing the mathematical strain experienced by `whfast` inside-out Jacobi coordinate hierarchies. | Extremely fast; allows timesteps of $10^4$ to $10^5$ years because orbital timescales are integrated out. |

### 2.1 High-Fidelity Initialization and Boundary-Validation Phase
Because a full secular shepherding integration requires $\sim 10\text{–}100$ Myr to manifest (the perturber's period at $550\text{ AU}$ is $\approx 12,900$ years, meaning 20,000 years spans only $\approx 1.55$ orbits), the 20,000-year integration length is not a steady-state secular configuration. Instead, it is explicitly defined as a **High-Fidelity Initialization and Boundary-Validation Phase** ($\Delta t \to t_{\text{max}}$ truncation check). This phase is mathematically used to confirm coordinate conversion stability under severe retrograde torque and ensure the code architecture is numerically converged before transferring it to high-performance computing (HPC) nodes.

---

## 3. Cartesian to Keplerian Orbital Element Transformations

To extract orbital parameters at each print step, we perform the following transformations relative to the central Sun (mass $M_\odot$):

1. **Specific Angular Momentum**:
   $$\mathbf{h} = \mathbf{r} \times \mathbf{v}, \quad h = |\mathbf{h}|$$

2. **Specific Orbital Energy**:
   $$\mathcal{E} = \frac{v^2}{2} - \frac{G M_\odot}{r}$_{}$

3. **Semi-Major Axis ($a$)**:
   $$a = -\frac{G M_\odot}{2 \mathcal{E}}$$

4. **Eccentricity Vector ($\mathbf{e}$)**:
   $$\mathbf{e} = \frac{1}{G M_\odot} \left[ \left(v^2 - \frac{G M_\odot}{r}\right)\mathbf{r} - (\mathbf{r} \cdot \mathbf{v})\mathbf{v} \right], \quad e = |\mathbf{e}|$$

5. **Inclination ($i$)**:
   $$i = \arccos\left(\frac{h_z}{h}\right)$$

6. **Longitude of Ascending Node ($\Omega$)**:
   Using the nodal vector $\mathbf{n} = \hat{\mathbf{z}} \times \mathbf{h} = (-h_y, h_x, 0)$:
   $$\Omega = \arctan2(h_x, -h_y)$$

7. **Argument of Perihelion ($\omega$)**:
   $$\omega = \arccos\left(\frac{\mathbf{n} \cdot \mathbf{e}}{n e}\right)$$
   If $e_z < 0$, then $\omega = 2\pi - \omega$. (For coplanar orbits, $\omega = \arctan2(e_y, e_x)$).

8. **Longitude of Perihelion ($\varpi$)**:
   $$\varpi = \Omega + \omega$$

---

## 4. Poincaré Cartesian Variables

To observe the secular perihelion alignment ($\Delta\varpi = \varpi - \varpi_9$) and nodal shepherding ($\Delta\Omega = \Omega - \Omega_9$) relative to the retrograde Planet Nine orbit, we project the elements into Cartesian-like Poincaré coordinates:

* **Eccentricity Poincaré coordinates**:
  $$h = e \sin(\Delta\varpi)$$
  $$k = e \cos(\Delta\varpi)$$

* **Inclination Poincaré coordinates**:
  $$p = \sin\left(\frac{i}{2}\right) \sin(\Delta\Omega)$$
  $$g = \sin\left(\frac{i}{2}\right) \cos(\Delta\Omega)$$

---

## 5. Atmosphere and Thermal Infrared Flux Equations

We model the direct thermal detectability of Planet Nine in the WISE W1 (3.4 $\mu$m) band, replicating the low-temperature atmosphere models of Fortney et al. (2016).

1. **Equilibrium Temperature ($T_{\text{eq}}$)**:
   $$T_{\text{eq}}(d) = T_\odot (1 - A)^{1/4} \sqrt{\frac{R_\odot}{2 d}} \approx 278.3\text{ K} \times (1 - A)^{1/4} \left(\frac{1\text{ AU}}{d}\right)^{1/2}$$
   Where $A = 0.3$ is the geometric albedo and $d$ is the heliocentric distance.

2. **Effective Temperature ($T_{\text{eff}}$)**:
   $$T_{\text{eff}}(d) = \left( T_{\text{int}}^4 + T_{\text{eq}}(d)^4 \right)^{1/4}$$
   Where $T_{\text{int}} \approx 40$ K represents the internal cooling temperature at 4.5 Gyr.

3. **WISE W1 (3.4 $\mu$m) Apparent Magnitude**:
   Under stratospheric methane freeze-out (depleted case, mixing ratio $\le 10^{-6}$), the atmospheric opacity decreases, increasing W1-band emission:
   $$m_{\text{W1}}(d, T_{\text{eff}}) = 17.5 + 5 \log_{10}\left(\frac{d}{700\text{ AU}}\right) - 15 \log_{10}\left(\frac{T_{\text{eff}}}{40\text{ K}}\right)$$
   For the standard methane-rich atmosphere:
   $$m_{\text{W1\_rich}}(d, T_{\text{eff}}) = 23.5 + 5 \log_{10}\left(\frac{d}{700\text{ AU}}\right) - 15 \log_{10}\left(\frac{T_{\text{eff}}}{40\text{ K}}\right)$$

4. **Infrared Flux Density ($F_{\mu\text{Jy}}$)**:
   $$F_{\mu\text{Jy}} = 10^6 \times F_0 \times 10^{-m_{\text{W1}} / 2.5}$$
   Where $F_0 = 309.54$ Jy is the zero-point flux of the WISE W1 band.
