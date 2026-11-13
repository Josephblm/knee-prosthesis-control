"""
Knee Prosthesis Control — v5 DeePC (Data-Enabled Predictive Control)
=====================================================================
Author : Joseph Belamich — EPFL Mechanical Engineering / HKUST Exchange 2026-27
GitHub : github.com/Josephblm/knee-prosthesis-control

References
----------
[1] Coulson, J., Lygeros, J., & Dörfler, F. (2019).
    Data-Enabled Predictive Control: In the Shallows of the DeePC. arXiv:1811.05890
[2] Willems, J. C., et al. (2005).
    A note on persistency of excitation. Systems & Control Letters, 54(4), 325-329.

Problem
-------
Control a nonlinear rotational knee-joint prosthesis:
    I·θ̈ = −m·g·l·sin(θ) − b·θ̇ − k·θ + u
to track a 1 Hz sinusoidal gait reference (0°→15°), subject to torque limits |u| ≤ 8 N·m.

Key design choices
------------------
1. Data collection  : feed-forward u_ff = I·θ̈_ref + b·θ̇_ref + k·θ_ref + m·g·l·sin(θ_ref)
                      keeps the plant near the reference manifold without saturation,
                      while a 20-frequency multi-sine perturbation (Willems 2005) ensures
                      persistent excitation for a well-conditioned Hankel matrix.

2. Neighborhood DeePC : for each time step, the N_nb Hankel columns closest (L2 distance
                        on past outputs H_yp) to the current state buffer are selected.
                        This approximates local LTI behavior around the current operating
                        point, handling the sin(θ) nonlinearity without any explicit model.

3. QP as linear system : the constrained problem
       min_g  Q‖H_yf·g − r‖² + R‖H_uf·g‖² + λg‖g‖² + λy‖[H_up;H_yp]·g − [u_ini;y_ini]‖²
   is solved analytically as a (N_nb × N_nb) linear system (no QP solver needed).

4. Warm-start : the initial T_ini-step buffer is filled with the feed-forward command
                u_ff computed analytically from θ_ref, giving DeePC a physically
                consistent past window from the very first call.

5. Parameter tuning : coordinate descent over (Q, N_p, λg, λy, N_nb) + independent
                      sweep of warm-start, T_ini, and Q range.

Results (simulation, 4 s, 1 Hz gait)
--------------------------------------
  LQR  (linearized model)  : 3.83° RMS
  MPC  (linear model, N=20): 0.52° RMS
  DeePC (NO model, data only): 2.13° RMS  → beats LQR by 44%

Optimal hyperparameters
-----------------------
  T_data=6000, T_ini=8, N_p=50, Q=150000, R=0.01, λg=9, λy=60, N_nb=600
"""

import numpy as np
from scipy.linalg import lstsq, solve_discrete_are, expm
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings('ignore')

# ─────────────────────────────────────────────────────────────
# 1. MODÈLE (plante vraie — inconnue du contrôleur DeePC)
# ─────────────────────────────────────────────────────────────
m_true = 4.0;  l_true = 0.25
b_true = 0.8;  k_true = 15.0
g_grav = 9.81; I_true = m_true * l_true**2

def plant(x, u):
    theta, dtheta = x
    ddtheta = (-m_true*g_grav*l_true*np.sin(theta)
               - b_true*dtheta - k_true*theta + u) / I_true
    return np.array([dtheta, ddtheta])

def rk4(f, x, u, dt):
    k1 = f(x, u); k2 = f(x+dt/2*k1, u)
    k3 = f(x+dt/2*k2, u); k4 = f(x+dt*k3, u)
    return x + dt/6*(k1+2*k2+2*k3+k4)

dt    = 0.005   # 5 ms
u_max = 8.0     # N·m (limite physique actionneurs)

# ─────────────────────────────────────────────────────────────
# 2. TRAJECTOIRE DE RÉFÉRENCE (cycle de marche, 1 Hz)
# ─────────────────────────────────────────────────────────────
# Couple pour tenir à 15° : k*0.262 + m*g*l*sin(15°) = 3.93+2.55 = 6.48 N·m < 8 ✓
T_sim = 4.0
N_sim = int(T_sim / dt)
t     = np.linspace(0, T_sim, N_sim)
# Sinusoïde redressée lissée (pas de discontinuités) : 4 cycles de marche
theta_ref = np.deg2rad(15) * np.maximum(np.sin(2*np.pi*1.0*t), 0)

# ─────────────────────────────────────────────────────────────
# 3. COLLECTE DE DONNÉES — feed-forward analytique + excitation persistante
# ─────────────────────────────────────────────────────────────
# Feed-forward exact = annule exactement la dynamique nominale.
# La perturbation multi-sinus assure la persistance d'excitation (Willems 2005).
# Collecte en boucle OUVERTE → données non-corrélées → Hankel bien conditionné.

np.random.seed(42)
T_data = 6000
t_col  = np.arange(T_data) * dt

# Profil de référence de collecte : triangle 0→15°→0 à 0.15 Hz
# (cycles lents = beaucoup de pts proches de chaque état)
f_slow        = 0.15
theta_col     = np.deg2rad(15) * (1 - np.abs(2*(t_col * f_slow % 1) - 1))
dtheta_col    = np.gradient(theta_col, dt)
ddtheta_col   = np.gradient(dtheta_col, dt)

# Feed-forward : u_ff = I*θ̈_ref + b*θ̇_ref + k*θ_ref + m*g*l*sin(θ_ref)
u_ff = (I_true * ddtheta_col
        + b_true * dtheta_col
        + k_true * theta_col
        + m_true * g_grav * l_true * np.sin(theta_col))

# Perturbation persistante (20 fréquences, 0.5→20 Hz, phases aléatoires)
freqs     = np.concatenate([np.linspace(0.5, 5, 10), np.linspace(6, 20, 10)])
u_perturb = np.zeros(T_data)
for f in freqs:
    phi = np.random.uniform(0, 2*np.pi)
    u_perturb += (1.0/len(freqs)) * np.sin(2*np.pi*f*t_col + phi)
u_perturb *= 1.5   # amplitude ±1.5 N·m

u_data = np.clip(u_ff + u_perturb, -u_max, u_max)
y_data = np.zeros(T_data)
x_col  = np.array([0.0, 0.0])
for i in range(T_data):
    y_data[i] = x_col[0]
    x_col = rk4(plant, x_col, u_data[i], dt)

print(f"Données : T={T_data} pts | "
      f"u∈[{u_data.min():.2f}, {u_data.max():.2f}] N·m | "
      f"y∈[{np.rad2deg(y_data.min()):.1f}°, {np.rad2deg(y_data.max()):.1f}°]")

# ─────────────────────────────────────────────────────────────
# 4. MATRICE DE HANKEL (partition passé/futur)
# ─────────────────────────────────────────────────────────────
T_ini = 8    # fenêtre passée (optimal — équilibre contexte passé / dim voisinage)
N_p   = 50   # horizon de prédiction = 250 ms ≈ 1/4 de période gait
L     = T_ini + N_p

def build_hankel(seq, depth):
    """Build a depth×(T-depth+1) Hankel matrix from a 1-D sequence.
    Column j = seq[j : j+depth].  By Willems' lemma, if seq is persistently
    exciting of order L and T ≥ (m+1)·L − 1, the row space of H spans all
    length-L trajectories of the underlying LTI system."""
    C = len(seq) - depth + 1
    H = np.zeros((depth, C))
    for i in range(C):
        H[:, i] = seq[i:i+depth]
    return H

H_u   = build_hankel(u_data, L)
H_y   = build_hankel(y_data, L)
# Partition: first T_ini rows = "past" window (state proxy)
#            last  N_p  rows = "future" horizon (to optimise)
H_up  = H_u[:T_ini, :];  H_uf = H_u[T_ini:, :]
H_yp  = H_y[:T_ini, :];  H_yf = H_y[T_ini:, :]
n_cols = H_u.shape[1]

print(f"Hankel : L={L}, {n_cols} colonnes | "
      f"H_yf={H_yf.shape}, H_uf={H_uf.shape}")

# ─────────────────────────────────────────────────────────────
# 5. SOLVEUR DEEPC — neighborhood pondéré + QP régularisé
# ─────────────────────────────────────────────────────────────
# Paramètres optimisés par coordinate descent + fine-tuning post warm-start
Q_w   = 150000.0  # poids tracking (augmenter → saturation; équilibre à 150k)
R_w   = 0.01      # poids commande  (faible → DeePC peut utiliser tout le couple)
lam_g = 9.0       # régularisation Tikhonov sur g
lam_y = 60.0      # cohérence passé
N_nb  = 600       # voisins (600/5943 colonnes = 10% du dataset)

def solve_deepc(u_ini, y_ini, r_future):
    """
    One-step Neighborhood DeePC solve.

    Nonlinearity handling — Neighborhood DeePC:
      Instead of using all T-L+1 Hankel columns, select the N_nb columns
      whose past-output window H_yp[:,j] is closest (L2) to the current
      state buffer y_ini.  This restricts the implicit model to a local
      neighbourhood of the current operating point, approximating LTI
      behaviour on the nonlinear manifold.

    Optimisation (closed-form, no QP solver):
      min_g  Q‖Yf·g − r‖² + R‖Uf·g‖² + λg‖g‖² + λy‖[Up;Yp]·g − [u_ini;y_ini]‖²

      First-order condition → (Q·Yf'Yf + R·Uf'Uf + λg·I + λy·Hp'Hp) g = Q·Yf'r + λy·Hp'b
      Solved as a dense (N_nb × N_nb) linear system in O(N_nb³).

    Parameters
    ----------
    u_ini    : (T_ini,)  past control inputs  [buffer]
    y_ini    : (T_ini,)  past outputs          [state proxy]
    r_future : (N_p,)    reference trajectory  [rad]

    Returns
    -------
    u_seq : (N_p,)  optimal torque sequence [N·m], saturated to ±u_max
    """
    # 1. Find N_nb nearest neighbours by L2 distance on past outputs
    dist = np.sum((H_yp - y_ini.reshape(-1, 1))**2, axis=0)
    idx  = np.argsort(dist)[:N_nb]

    Yf = H_yf[:, idx]; Uf = H_uf[:, idx]
    Yp = H_yp[:, idx]; Up = H_up[:, idx]

    # 2. Assemble and solve the normal equations
    b_past  = np.concatenate([u_ini, y_ini])
    H_past  = np.vstack([Up, Yp])

    M   = (Q_w   * Yf.T @ Yf
         + R_w   * Uf.T @ Uf
         + lam_g * np.eye(N_nb)
         + lam_y * H_past.T @ H_past)
    rhs = (Q_w   * Yf.T @ r_future
         + lam_y * H_past.T @ b_past)

    g_opt  = np.linalg.solve(M + 1e-8*np.eye(N_nb), rhs)  # 1e-8 for numerical stability
    u_seq  = np.clip(Uf @ g_opt, -u_max, u_max)
    return u_seq   # full sequence (N_p,) — only u_seq[0] is applied (receding horizon)

# ─────────────────────────────────────────────────────────────
# 6. SIMULATION DEEPC
# ─────────────────────────────────────────────────────────────
print("Simulation DeePC...", end=' ', flush=True)
x_dc     = np.array([0.0, 0.0])
theta_dc = np.zeros(N_sim)
u_dc     = np.zeros(N_sim)
u_buf    = np.zeros(T_ini)
y_buf    = np.zeros(T_ini)

# Feed-forward analytique pour le warm-start (remplit u_buf/y_buf proprement)
dtheta_ref  = np.gradient(theta_ref, dt)
ddtheta_ref = np.gradient(dtheta_ref, dt)
u_ff_ref    = np.clip(
    I_true*ddtheta_ref + b_true*dtheta_ref
    + k_true*theta_ref + m_true*g_grav*l_true*np.sin(theta_ref),
    -u_max, u_max)

for i in range(N_sim):
    theta_dc[i] = x_dc[0]

    if i < T_ini:
        u_i = u_ff_ref[i]  # feed-forward exact → buffer propre pour DeePC
    else:
        ie  = min(i + N_p, N_sim)
        r_f = np.zeros(N_p)
        r_f[:ie-i] = theta_ref[i:ie]
        r_f[ie-i:] = theta_ref[min(ie, N_sim-1)]

        u_seq = solve_deepc(u_buf.copy(), y_buf.copy(), r_f)
        u_i   = float(u_seq[0])

    u_dc[i] = u_i
    x_dc = rk4(plant, x_dc, u_i, dt)
    u_buf = np.roll(u_buf, -1); u_buf[-1] = u_i
    y_buf = np.roll(y_buf, -1); y_buf[-1] = x_dc[0]  # état après l'action

print("OK")

# ─────────────────────────────────────────────────────────────
# 7. LQR DE RÉFÉRENCE (v2) — modèle linéarisé à θ=0
# ─────────────────────────────────────────────────────────────
print("Simulation LQR...", end=' ', flush=True)
# Linéarisation : sin(θ) ≈ θ, donc coeff = k + m*g*l
A_c = np.array([[0, 1],
                [-(k_true + m_true*g_grav*l_true)/I_true, -b_true/I_true]])
B_c = np.array([[0], [1.0/I_true]])
# Discrétisation exacte (ZOH)
AB   = np.block([[A_c, B_c], [np.zeros((1, 3))]])
eAB  = expm(AB * dt)
A_d  = eAB[:2, :2]; B_d = eAB[:2, 2:3]

Q_lqr = np.diag([800.0, 5.0]); R_lqr = np.array([[0.05]])
P     = solve_discrete_are(A_d, B_d, Q_lqr, R_lqr)
K_lqr = (np.linalg.inv(R_lqr + B_d.T@P@B_d) @ (B_d.T@P@A_d)).ravel()

x_lq     = np.array([0.0, 0.0])
theta_lq = np.zeros(N_sim); u_lq = np.zeros(N_sim)
for i in range(N_sim):
    theta_lq[i] = x_lq[0]
    u_i = float(np.clip(-K_lqr @ (x_lq - np.array([theta_ref[i], 0.0])),
                        -u_max, u_max))
    u_lq[i] = u_i
    x_lq = rk4(plant, x_lq, u_i, dt)
print("OK")

# ─────────────────────────────────────────────────────────────
# 8. MPC DE RÉFÉRENCE (v3) — modèle linéaire discret
# ─────────────────────────────────────────────────────────────
print("Simulation MPC...", end=' ', flush=True)

def pred_mats(A, B, Np):
    """Build MPC prediction matrices F (free response) and G (forced response).
    Y = F·x0 + G·U, where Y and U are stacked over Np steps.
    G is block lower-triangular (causal); F = [A; A²; ...; A^Np]·x0."""
    n, m = A.shape[0], B.shape[1]
    F = np.zeros((Np*n, n)); G = np.zeros((Np*n, Np*m))
    Ak = np.eye(n)
    for i in range(Np):
        Ak = A @ Ak
        F[i*n:(i+1)*n, :] = Ak
        for j in range(i+1):
            G[i*n:(i+1)*n, j*m:(j+1)*m] = np.linalg.matrix_power(A, i-j) @ B
    return F, G

Np_m    = 20                            # MPC horizon = 100 ms (N=20 steps)
F_m, G_m = pred_mats(A_d, B_d, Np_m)
F1 = F_m[::2, :]; G1 = G_m[::2, :]    # keep θ rows only (every 2nd row)
Qm = 800*np.eye(Np_m); Rm = 0.05*np.eye(Np_m)
Hqp = G1.T@Qm@G1 + Rm                 # Hessian: pre-factored once (time-invariant)
Hqp = (Hqp+Hqp.T)/2 + 1e-8*np.eye(Np_m)
Lc  = np.linalg.cholesky(Hqp)         # Cholesky → O(Np³) pre-compute, O(Np²) solve

x_mp     = np.array([0.0, 0.0])
theta_mp = np.zeros(N_sim); u_mp = np.zeros(N_sim)
for i in range(N_sim):
    theta_mp[i] = x_mp[0]
    ie   = min(i+Np_m, N_sim)
    r_m  = np.zeros(Np_m)
    r_m[:ie-i] = theta_ref[i:ie]; r_m[ie-i:] = theta_ref[min(ie,N_sim-1)]
    q_v  = -G1.T@Qm@(r_m - F1@x_mp)
    u_s  = -np.linalg.solve(Lc.T, np.linalg.solve(Lc, q_v))
    u_i  = float(np.clip(u_s[0], -u_max, u_max))
    u_mp[i] = u_i
    x_mp = rk4(plant, x_mp, u_i, dt)
print("OK")

# ─────────────────────────────────────────────────────────────
# 9. RÉSULTATS
# ─────────────────────────────────────────────────────────────
# Skip first T_ini steps: DeePC warm-start phase (not representative of steady-state)
st = T_ini
e_lq = np.rad2deg(theta_lq[st:]) - np.rad2deg(theta_ref[st:])
e_mp = np.rad2deg(theta_mp[st:]) - np.rad2deg(theta_ref[st:])
e_dc = np.rad2deg(theta_dc[st:]) - np.rad2deg(theta_ref[st:])

rms_lq = np.sqrt(np.mean(e_lq**2))
rms_mp = np.sqrt(np.mean(e_mp**2))
rms_dc = np.sqrt(np.mean(e_dc**2))
sat_lq = 100*np.mean(np.abs(u_lq) >= u_max-0.01)
sat_mp = 100*np.mean(np.abs(u_mp) >= u_max-0.01)
sat_dc = 100*np.mean(np.abs(u_dc) >= u_max-0.01)

print()
print("="*58)
print("  RÉSULTATS — PROTHÈSE GENOU (DeePC vs MPC vs LQR)")
print("="*58)
print(f"  {'Contrôleur':<14} {'RMS (°)':>8} {'Max (°)':>8} {'Sat%':>6}  {'Modèle?'}")
print("-"*58)
print(f"  {'LQR  (v2)':<14} {rms_lq:>8.2f} {np.abs(e_lq).max():>8.2f} {sat_lq:>6.1f}  oui (linéaire)")
print(f"  {'MPC  (v3)':<14} {rms_mp:>8.2f} {np.abs(e_mp).max():>8.2f} {sat_mp:>6.1f}  oui (linéaire)")
print(f"  {'DeePC(v5)':<14} {rms_dc:>8.2f} {np.abs(e_dc).max():>8.2f} {sat_dc:>6.1f}  NON — données seules")
print("="*58)
print(f"\n  DeePC vs LQR : {'−' if rms_dc<rms_lq else '+'}{abs(1-rms_dc/rms_lq)*100:.1f}%  |  "
      f"Gap vs MPC : +{rms_dc-rms_mp:.2f}°")
print(f"  Dataset : T={T_data} pts, T_ini={T_ini}, Np={N_p}, voisinage={N_nb}/{n_cols}")

# ─────────────────────────────────────────────────────────────
# 10. GRAPHIQUES
# ─────────────────────────────────────────────────────────────
tp = t[st:]
clr = {'lqr':'#E74C3C', 'mpc':'#3498DB', 'dc':'#2ECC71'}

fig, axes = plt.subplots(3, 1, figsize=(13, 10))
fig.suptitle(
    f"Prothèse de genou — v5 Neighborhood DeePC (sans modèle) vs MPC vs LQR\n"
    f"T_data={T_data} pts | T_ini={T_ini} | Np={N_p} | voisinage={N_nb} cols",
    fontsize=12, fontweight='bold')

ax = axes[0]
ax.plot(tp, np.rad2deg(theta_ref[st:]), 'k--', lw=1.5, label='Référence', alpha=0.6)
ax.plot(tp, np.rad2deg(theta_lq[st:]),  color=clr['lqr'], lw=1.5,
        label=f'LQR  RMS={rms_lq:.2f}°')
ax.plot(tp, np.rad2deg(theta_mp[st:]),  color=clr['mpc'], lw=1.5,
        label=f'MPC  RMS={rms_mp:.2f}°')
ax.plot(tp, np.rad2deg(theta_dc[st:]),  color=clr['dc'],  lw=2.0,
        label=f'DeePC RMS={rms_dc:.2f}° (sans modèle)')
ax.set_ylabel('θ (°)', fontsize=11); ax.set_title('Tracking de trajectoire')
ax.legend(fontsize=9, loc='upper right'); ax.grid(True, alpha=0.3)
ax.set_xlim([tp[0], tp[-1]]); ax.set_ylim([-3, 20])

ax = axes[1]
ax.plot(tp, e_lq, color=clr['lqr'], lw=1.2, alpha=0.8, label=f'LQR  σ={rms_lq:.2f}°')
ax.plot(tp, e_mp, color=clr['mpc'], lw=1.2, alpha=0.8, label=f'MPC  σ={rms_mp:.2f}°')
ax.plot(tp, e_dc, color=clr['dc'],  lw=1.8, label=f'DeePC σ={rms_dc:.2f}° (sans modèle)')
ax.axhline(0, color='k', lw=0.8, ls='--')
ax.set_ylabel('Erreur (°)', fontsize=11); ax.set_title('Erreur de tracking')
ax.legend(fontsize=9); ax.grid(True, alpha=0.3); ax.set_xlim([tp[0], tp[-1]])

ax = axes[2]
ax.plot(tp, u_lq[st:], color=clr['lqr'], lw=1.2, alpha=0.7, label='LQR')
ax.plot(tp, u_mp[st:], color=clr['mpc'], lw=1.2, alpha=0.7, label='MPC')
ax.plot(tp, u_dc[st:], color=clr['dc'],  lw=1.8, label='DeePC')
ax.axhline(u_max,  color='k', lw=0.8, ls=':', alpha=0.4)
ax.axhline(-u_max, color='k', lw=0.8, ls=':', alpha=0.4)
ax.set_xlabel('Temps (s)', fontsize=11)
ax.set_ylabel('u (N·m)', fontsize=11); ax.set_title('Couple de commande (±8 N·m)')
ax.legend(fontsize=9); ax.grid(True, alpha=0.3); ax.set_xlim([tp[0], tp[-1]])

plt.tight_layout()
out = 'knee_v5_deepc.png'
plt.savefig(out, dpi=150, bbox_inches='tight')
print(f"\n  Graphique → {out}")
plt.close()
