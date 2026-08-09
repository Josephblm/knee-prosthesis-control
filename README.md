# knee-prosthesis-control
Simulation of an LQR controller for a variable-stiffness knee prosthesis in MATLAB
# Bionic Knee Prosthesis — LQR Control Simulation

**Joseph Belamich | EPFL Mechanical Engineering | Exchange Student at HKUST**

## Overview
This project simulates an LQR-based optimal controller for a bionic knee 
prosthesis, inspired by recent work on variable-stiffness prosthetic devices 
(Herr et al., 2024). The goal is to reproduce a natural gait trajectory while 
maintaining robustness to sensor noise and adapting joint stiffness to the 
phase of walking.

## Motivation
Millions of people worldwide live with lower-limb amputations. Current 
prosthetic knees struggle to reproduce the natural biomechanics of human 
walking. This simulation explores how optimal control theory can help bridge 
that gap — a problem I want to pursue at the graduate level.

## System Model
The knee joint is modeled as a rotational mass-spring-damper system:

**I·θ̈ + b·θ̇ + k·θ = τ**

Where:
- I = 0.5 kg·m² — moment of inertia of the lower leg
- b = 2.0 N·m·s/rad — joint damping
- k = variable N·m/rad — joint stiffness (phase-dependent)
- τ = motor torque command (N·m)

The system is expressed in state-space form with state vector x = [θ, θ̇].

## Controller Design
An LQR (Linear Quadratic Regulator) controller is designed to minimize:

**J = ∫ (xᵀQx + uᵀRu) dt**

With Q = diag([2000, 50]) and R = 0.01, prioritizing tracking accuracy 
over energy consumption.

## Results

### 1. LQR Tracking — Single Gait Cycle
![LQR Control](fig2_lqr_control.png)
The prosthesis accurately tracks the natural gait reference trajectory.

### 2. Comparison: With vs Without Controller
![Comparison](fig3_comparison.png)
Without control, the knee remains at 0°. The LQR drives it to follow 
the natural walking pattern.

### 3. Multiple Gait Cycles
![Multiple Cycles](fig4_multiple_cycles.png)
The controller maintains stable tracking over 4 consecutive steps.

### 4. Robustness to Sensor Noise
![Noise Robustness](fig5_noise_robustness.png)
Performance is maintained under realistic sensor noise (σ = 0.5°).

### 5. LQR vs PID Comparison
![LQR vs PID](fig5_comparison_pid.png)
The LQR outperforms a classical PID controller, with significantly 
reduced phase lag.

### 6. Variable Stiffness Control
![Variable Stiffness](fig6_variable_stiffness.png)
Inspired by Herr et al. (2024): stiffness is modulated between stance 
(k=20 N·m/rad) and swing (k=3 N·m/rad) phases, reducing motor effort 
while maintaining tracking accuracy.

## Key Findings
- LQR achieves near-perfect trajectory tracking with <2° RMS error
- Variable stiffness reduces peak motor torque by ~40% vs fixed stiffness
- Controller remains stable under sensor noise up to σ = 2°

## References
- Herr, H. et al. (2024). Variable-stiffness prosthesis improves 
  biomechanics of walking across speeds. Scientific Reports.
- Hogan, N. (1985). Impedance Control: An Approach to Manipulation. 
  Journal of Dynamic Systems.

## Tools
MATLAB R2024 — Control System Toolbox (lqr, ode45)

## Author
Joseph Belamich — joseph.belamich@epfl.ch  
EPFL Mechanical Engineering, GPA 5.5/6.0  
Exchange Student, HKUST 2026-2027
