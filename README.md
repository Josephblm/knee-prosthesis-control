# Bionic Knee Prosthesis — Control Simulation Portfolio

**Joseph Belamich | EPFL Mechanical Engineering | Exchange at HKUST 2026–2027**

---

## Overview

Progressive simulation of a bionic knee prosthesis controller, evolving from classical optimal control (LQR) to data-driven model-free control (DeePC). The system models a knee joint as a nonlinear rotational mass-spring-damper with actuator torque limits.

The project is structured as a learning trajectory: each version adds a layer of complexity and realism.

---

## System Model

The knee joint dynamics:

**I·θ̈ = −m·g·l·sin(θ) − b·θ̇ − k·θ + u**

| Parameter | Value | Description |
|---|---|---|
| m | 4.0 kg | Shank mass |
| l | 0.25 m | CoM distance |
| I = ml² | 0.25 kg·m² | Moment of inertia |
| b | 0.8 N·m·s/rad | Joint damping |
| k | 15.0 N·m/rad | Joint stiffness |
| u_max | ±8 N·m | Torque saturation |

Reference trajectory: 15° sinusoidal gait at 1 Hz (physically achievable with u_max given the nonlinear gravity term).

---

## Controllers — Comparative Results

| Version | Controller | RMS Error | Max Error | Model required? |
|---|---|---|---|---|
| v2 | **LQR** (linearized) | 3.83° | 7.71° | Yes (linear) |
| v3 | **MPC** (receding horizon) | 0.52° | 2.54° | Yes (linear) |
| **v5** | **DeePC** (data-driven) | **2.13°** | **7.75°** | **NO — data only** |

DeePC beats LQR by **44%** using no physical model whatsoever — only 6000 input-output data samples.

---

## v5 — DeePC (Data-Enabled Predictive Control)

Based on [Coulson et al., 2019](https://arxiv.org/abs/1811.05890) and the Willems fundamental lemma.

### Key idea
Instead of identifying a parametric model, DeePC builds a **Hankel matrix** from open-loop input-output data. At each timestep, it finds a linear combination of data trajectories that matches the current system state ("neighborhood"), then optimizes future inputs over this local linear approximation.

This handles nonlinearity via **Neighborhood DeePC**: only the K nearest Hankel columns (by L2 distance on past outputs) are used, approximating local LTI behavior.

### Data collection strategy
- **Feed-forward analytical**: u_ff = I·θ̈_ref + b·θ̇_ref + k·θ_ref + m·g·l·sin(θ_ref) — covers 0°→15° without saturation
- **Multi-sinus perturbation** (20 frequencies, 0.5–20 Hz, ±1.5 N·m) — ensures persistent excitation (Willems 2005)
- Open-loop collection → uncorrelated data → well-conditioned Hankel

### Optimization problem (solved as linear system)

```
min_g  Q·‖H_yf·g − r‖² + R·‖H_uf·g‖² + λg·‖g‖² + λy·‖[H_up; H_yp]·g − [u_ini; y_ini]‖²
```

### Optimal parameters (coordinate descent + full sweep, converged)

| Parameter | Value | Role |
|---|---|---|
| T_data | 6000 samples | Dataset size |
| T_ini | 8 steps | Past window (state proxy) |
| N_p | 50 steps | Prediction horizon (250 ms lookahead) |
| Q_w | 150000 | Tracking weight |
| R_w | 0.01 | Control effort weight |
| λg | 9 | Tikhonov regularization |
| λy | 60 | Past consistency |
| N_nb | 600 | Neighborhood size (out of 5943 columns) |
| warm-start | feed-forward | Fills initial buffer with u_ff instead of zeros |

---

## Files

| File | Description |
|---|---|
| `knee_prosthesis.m` | MATLAB: LQR + variable stiffness (v1–v2) |
| `knee_v5_deepc.py` | Python: LQR + MPC + DeePC full comparison (v5) |
| `knee_v5_deepc.png` | Result figure: tracking, error, torque |

---

## References

- Coulson, J., Lygeros, J., & Dörfler, F. (2019). *Data-Enabled Predictive Control: In the Shallows of the DeePC*. arXiv:1811.05890
- Willems, J. C., et al. (2005). *A note on persistency of excitation*. Systems & Control Letters.
- Herr, H., et al. (2024). Variable-stiffness prosthetic knee. MIT Media Lab.
- Hogan, N. (1985). *Impedance control: An approach to manipulation*. ASME J. Dyn. Sys.
