"""
Figure 4 – Schematic mapping from operator spectrum to critical‑line zeros
===========================================================================
Left: simulated eigenvalues λ_n^2 of the squared operator Δ_K.  The
horizontal dashed line marks the threshold 1/4.  Right: once the identity
s(1‑s) = λ_n^2 is imposed, each pair (λ_n, ‑λ_n) is mapped onto a
conjugate pair of zeros on Re(s) = 1/2.  Output: fig4_zero_mapping.pdf.
"""
import numpy as np
import matplotlib.pyplot as plt

# ---------- Simulate eigenvalues λ_n^2 ----------
n = np.arange(1, 21)
# asymptotic law: λ_n ~ 2π n / ln n
lambda_sq = (2 * np.pi * n / np.log(np.maximum(n, 2)))**2

# ---------- Map to the s‑plane ----------
valid = lambda_sq > 0.25
t = np.sqrt(lambda_sq[valid] - 0.25)
s_upper = 0.5 + 1j * t
s_lower = 0.5 - 1j * t

# ---------- Figure ----------
plt.rcParams.update({'font.size': 11, 'font.family': 'serif'})
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

# Left panel: eigenvalues
ax1.plot(n, lambda_sq, 'bo', markersize=5, label=r'$\lambda_n^2$')
ax1.axhline(y=0.25, color='gray', linestyle='--', linewidth=1,
            label='Threshold 1/4')
ax1.set_xlabel('n')
ax1.set_ylabel(r'$\lambda_n^2$')
ax1.set_title('Spectrum of ' + r'$\Delta_K$')
ax1.legend(frameon=False)
ax1.grid(True, alpha=0.3)

# Right panel: critical‑line zeros
ax2.plot(np.real(s_upper), np.imag(s_upper), 'r.', markersize=5,
         label='Non‑trivial zeros')
ax2.plot(np.real(s_lower), np.imag(s_lower), 'r.', markersize=5)
ax2.axvline(x=0.5, color='gray', linestyle='--', linewidth=1,
            label='Critical line')
ax2.set_xlabel('Re(s)')
ax2.set_ylabel('Im(s)')
ax2.set_title('Mapping to zeta zeros')
ax2.legend(frameon=False)
ax2.set_xlim([0, 1])
ax2.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('fig4_zero_mapping.pdf', dpi=300, bbox_inches='tight')
plt.show()