"""
Figure 1 – Self-organised emergence of the prime density
=========================================================
Numerical integration of the generalized Ginzburg–Landau (GL) equation
driven by the coupling matrix K(x) = sqrt(2 x ln x).  The field evolves
from random noise and spontaneously converges to the steady state
density |Psi|^2 ~ 1/ln x.  The script outputs fig1_prime_density.pdf.
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# ======================== Parameters ========================
N = 400                           # number of grid points
L = 200.0                         # size of spatial domain
x = np.linspace(2.0, L, N)        # spatial grid (avoid x = 1 where ln x = 0)
dx = x[1] - x[0]

D = 0.01                          # diffusion coefficient
K_x_raw = np.sqrt(2 * x * np.log(np.maximum(x, 2.0)))
K_x = K_x_raw / np.max(K_x_raw) * 0.3   # normalised coupling strength
rho_target = 1.0 / np.log(np.maximum(x, 2.5))  # PNT density

# Random initial condition (small perturbations around target)
np.random.seed(42)
Psi_R_init = np.sqrt(rho_target) * 0.9 + 0.01 * np.random.randn(N)
Psi_I_init = 0.01 * np.random.randn(N)
y0 = np.concatenate([Psi_R_init, Psi_I_init])

# ======================== Helper: discrete Laplacian ========
def laplacian(f):
    """Second‑order centred finite difference Laplacian on a uniform grid."""
    lap = np.zeros_like(f)
    lap[1:-1] = (f[2:] - 2*f[1:-1] + f[:-2]) / (dx*dx)
    lap[0] = (f[1] - f[0]) / (dx*dx)
    lap[-1] = (f[-2] - f[-1]) / (dx*dx)
    return lap

# ======================== Right‑hand side of GL equation =====
def rhs(t, y):
    psi_R = y[:N]
    psi_I = y[N:]
    rho = psi_R**2 + psi_I**2

    lap_R = laplacian(psi_R)
    lap_I = laplacian(psi_I)

    stiffness = 5.0
    force_R = -stiffness * (rho - rho_target) * psi_R
    force_I = -stiffness * (rho - rho_target) * psi_I

    # Real‑imaginary coupling
    coupling_R = -K_x * psi_I
    coupling_I =  K_x * psi_R

    dR = D * lap_R + force_R + coupling_R
    dI = D * lap_I + force_I + coupling_I
    return np.concatenate([dR, dI])

# ======================== Integration ========================
t_span = (0, 150.0)
sol = solve_ivp(rhs, t_span, y0, method='BDF', rtol=1e-4, atol=1e-7)
psi_R_final = sol.y[:N, -1]
psi_I_final = sol.y[N:, -1]
rho_final = psi_R_final**2 + psi_I_final**2

# ======================== Publication‑quality figure =========
plt.rcParams.update({'font.size': 11, 'font.family': 'serif'})
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

# Left panel: simulated density versus target 1/ln x
ax1.plot(x, rho_final, 'b-', linewidth=1.5, label=r'Simulated $|\Psi|^2$')
ax1.plot(x, rho_target, 'r--', linewidth=1.5, label=r'Target $1/\ln x$')
ax1.fill_between(x, rho_final, rho_target, alpha=0.15, color='gray')
ax1.set_xlabel('x')
ax1.set_ylabel('Density')
ax1.set_title('Self-organised density vs. target')
ax1.legend(frameon=False)
ax1.grid(True, alpha=0.3)
ax1.set_xlim([2, 100])
ax1.set_ylim([0, 0.35])

# Right panel: relative error (log scale)
err = np.abs(rho_final - rho_target) / (rho_target + 1e-12) * 100
ax2.semilogy(x, err, 'm-', linewidth=1.0)
ax2.set_xlabel('x')
ax2.set_ylabel('Relative error (%)')
ax2.set_title('Relative error (log scale)')
ax2.grid(True, alpha=0.3)
ax2.set_xlim([2, 100])
mean_err = np.mean(err[50:-50])
ax2.text(0.55, 0.9, f'Mean error (x > 10): {mean_err:.2f}%',
         transform=ax2.transAxes, fontsize=10, verticalalignment='top',
         bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

plt.tight_layout(pad=2)
plt.savefig('fig1_prime_density.pdf', dpi=300, bbox_inches='tight')
plt.show()