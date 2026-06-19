# retrograde-p9

**A reproducible N-body pipeline and numerical convergence framework for retrograde Planet Nine configurations.**

This repository contains the integration pipeline, photometry parameterization, and convergence-boundary run described in the accompanying paper:

> *A Reproducible N-Body Pipeline and Numerical Convergence Framework for Retrograde Planet Nine Configurations* — Mayone Maha Rajan (Maha Strategies), Revision 3, June 2026. [Zenodo DOI to be added on deposit.]

---

## What this is, and what it is not

This is **infrastructure, not a result.** It does *not* claim that a retrograde perturber shepherds extreme trans-Neptunian objects (eTNOs). It provides a fully resolved active-point-mass integration framework and an initialization / boundary-validation run, designed to be scaled to high-performance computing (HPC).

The included 20,000-year run is a **convergence and stability check** — explicitly four to five orders of magnitude shorter than the 10–100 Myr secular shepherding timescale. Long-term shepherding cannot manifest within this run by construction. Any physical shepherding conclusion requires the HPC-scale integration this pipeline is built to enable.

If you are looking for a determination of shepherding efficiency, this repository does not contain one, and the paper says so plainly.

---

## Configuration studied

| Parameter | Value |
| :--- | :--- |
| Perturber mass `m₉` | 6.2 M⊕ |
| Semi-major axis `a₉` | 550 AU |
| Eccentricity `e₉` | 0.3 |
| Inclination `i₉` | 145° (retrograde) |
| Test particles `N` | 2,000 |
| Integration span | 20,000 years (≈ 1.55 perturber revolutions) |
| Integrator | REBOUND `ias15` (adaptive) |
| Random seed | 42 (pinned) |

Jupiter, Saturn, Uranus, and Neptune are modeled as **active moving point masses**, not secular "wires" — the pipeline resolves mean-motion resonances and close-encounter scattering rather than smoothing them away.

---

## ⚠️ Integrator: `ias15` only — the Leapfrog path is fenced off

**All reported results come from the REBOUND `ias15` adaptive integrator.**

The repository also contains a fixed-step vectorized NumPy **Leapfrog routine**. This exists **only as a dependency-failure fallback** so the pipeline can run when REBOUND is unavailable. It produced **none** of the results in the paper.

Leapfrog is a fixed-step symplectic method and **cannot faithfully resolve close-encounter scattering** (Neptune flybys) or mean-motion resonance capture. **Any run that uses the fallback path invalidates the resonance and flyby claims and must be flagged as such.** Do not present Leapfrog-path output as equivalent to `ias15` output.

---

## Requirements

- Python 3.x
- [REBOUND](https://rebound.readthedocs.io/) (provides the `ias15` integrator)
- NumPy, SciPy, Matplotlib

Pinned versions are in `requirements.txt` and the environment lockfile. **The lockfile is authoritative** — for exact reproduction, install from it rather than from loose version ranges.

```bash
# clone
git clone https://github.com/Maha-Strategies/retrograde-p9.git
cd retrograde-p9

# (recommended) create an isolated environment
python -m venv .venv
source .venv/bin/activate        # Windows: .venv\Scripts\activate

# install pinned dependencies
pip install -r requirements.txt
```

---

## Running the pipeline

```bash
python retrograde_p9_sim.py
```

With seed 42 and the pinned dependencies, the integration reproduces the paper's trajectories exactly. The run emits the phase-space and alignment figures:

- `poincare_eccentricity.png`
- `poincare_inclination.png`
- `perihelion_alignment_hist.png`

and the verified scalar metrics (orbital period `T₉ ≈ 12,899 yr`, revolution count `≈ 1.55`, W1 magnitude `m_W1 ≈ 16.98`, equilibrium temperature `≈ 11 K`, initial anti-aligned fraction `≈ 0.50`).

---

## ⚠️ Ephemerides: hardcoded for offline reproduction — re-pull before relying on them

The integration in the paper was produced in a **network-isolated sandbox**, so barycentric planet elements were **hardcoded at the module level** for offline stability. This is disclosed as a boundary condition, not hidden.

**If you are reproducing on a networked machine, re-pull live ephemerides from NASA JPL HORIZONS and confirm the hardcoded values before treating any element as authoritative.** See [`EPHEMERIDES.md`](./EPHEMERIDES.md) for exact instructions.

---

## Photometry

The infrared detectability parameterization (WISE W1, 3.4 μm) is anchored to the stratospheric methane freeze-out model of Fortney et al. (2016), *ApJL* 824, L25 (doi:10.3847/2041-8205/824/2/L25). It is a parameterization, **not** a line-by-line radiative-transfer solution, and should not be used as one.

---

## Provenance and verification

This repository follows a zero-fabrication, falsification-first discipline. Claims in the paper carry explicit tags:

- **[VERIFIED]** — internal arithmetic recomputed, or a cited parameterization confirmed against its source.
- **[SOURCED]** — drawn from cited literature in good faith; identifier not independently re-resolved.
- **[BOUNDARY]** — an explicit limit of the run, disclosed rather than worked around.

The three documented boundaries are: (1) timescale truncation (20,000 yr vs. the 10–100 Myr shepherding timescale); (2) offline hardcoded ephemerides; (3) parameterized photometry in place of radiative transfer. These are **conditions of the run, not defects** — they are the precise list of what would have to change before any physical claim about retrograde shepherding could be made.

---

## How to cite

If you use this pipeline, please cite the Zenodo deposit (DOI to be added) and this repository. Citation metadata is in [`CITATION.cff`](./CITATION.cff) once the DOI is minted.

## License

[Choose and add — e.g. MIT for code. The paper is deposited separately under CC-BY-4.0.]

## Acknowledgment of method

This pipeline was produced through human-directed AI synthesis: the architecture and audit protocol were designed and are owned by the author; an agentic AI instrument built the pipeline and drafted the manuscript under that direction; a second AI instrument performed editorial and provenance verification. Responsibility for all claims rests with the author.
