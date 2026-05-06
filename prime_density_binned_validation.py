"""
Prime density generator – binned π(x) validation
==================================================
The interval [2, 100] is divided into 10 equal sub‑intervals.  The
model‑predicted prime count (obtained by integrating the GL‑evolved
density |Ψ|^2) is compared with the exact count from sympy.primepi.
The script also reports the mean relative density error over the
interior of the domain.
"""
import numpy as np
from scipy.integrate import solve_ivp, cumulative_trapezoid
import sympy

# ======================== Parameters ========================
x_start, x_end = 2.0, 100.0
N = 500                          # grid points
t_end = 100.0                    # evolution time
D, kappa = 0.01, 5.0             # diffusion and feedback stiffness
seed = 42
bins = 10                        # number of sub‑intervals

# ======================== Grid and coupling ========================
x = np.linspace(x_start, x_end, N)
dx = x[1] - x[0]
K_raw = np.sqrt(2.0 * x * np.log(np.maximum(x, 2.0)))
K = K_raw / np.max(K_raw) * 0.3
rho_target = 1.0 / np.log(np.maximum(x, 2.5))

# ======================== Initial condition ========================
np.random.seed(seed)
psi_R = np.sqrt(rho_target) * 0.9 + 0.01 * np.random.randn(N)
psi_I = 0.01 * np.random.randn(N)
y0 = np.concatenate([psi_R, psi_I])

def laplacian(f):
    lap = np.zeros_like(f)
    lap[1:-1] = (f[2:] - 2.0 * f[1:-1] + f[:-2]) / (dx * dx)
    lap[0]    = (f[1] - f[0]) / (dx * dx)
    lap[-1]   = (f[-2] - f[-1]) / (dx * dx)
    return lap

def rhs(t, y):
    psi_R, psi_I = y[:N], y[N:]
    rho = psi_R**2 + psi_I**2
    lap_R, lap_I = laplacian(psi_R), laplacian(psi_I)
    force = -kappa * (rho - rho_target)
    dR = D * lap_R + force * psi_R - K * psi_I
    dI = D * lap_I + force * psi_I + K * psi_R
    return np.concatenate([dR, dI])

# ======================== Evolution ========================
sol = solve_ivp(rhs, (0, t_end), y0, method='BDF', rtol=1e-4, atol=1e-7)
psi_R_f, psi_I_f = sol.y[:N, -1], sol.y[N:, -1]
rho_final = psi_R_f**2 + psi_I_f**2

# ======================== Model π(x) ========================
integral = cumulative_trapezoid(rho_final, x, initial=0)
model_pi = 1.0 + integral            # π(2) = 1

def model_pi_at(x_val):
    return np.interp(x_val, x, model_pi)

# ======================== Binned validation ========================
bin_edges = np.linspace(x_start, x_end, bins + 1)
print(f"Binned π(x) comparison ({bins} intervals):")
print(f"{'Interval':>12} | {'True π':>6} | {'Model π':>8} | {'Abs. diff':>8}")
print("-" * 50)
for i in range(bins):
    a, b = bin_edges[i], bin_edges[i+1]
    true_count = int(sympy.primepi(float(b))) - int(sympy.primepi(float(a)))
    model_count = model_pi_at(float(b)) - model_pi_at(float(a))
    abs_diff = abs(true_count - model_count)
    print(f"[{a:4.1f}, {b:5.1f}] | {true_count:6d} | {model_count:8.2f} | "
          f"{abs_diff:8.2f}")

# ---------- Overall density error ----------
interior = slice(20, -20)
rel_err = np.mean(
    np.abs(rho_final[interior] - rho_target[interior])
    / rho_target[interior]
)
print(f"\nMean relative density error (interior): {rel_err*100:.3f}%")