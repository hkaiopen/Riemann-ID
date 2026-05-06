"""
Figure 2 – Coupling matrix K(x) and the prime density
=======================================================
The dynamic coupling matrix K(x) = sqrt(2 x ln x) (left axis, blue)
is plotted together with the PNT density 1/ln x (right axis, red).
The anti‑correlation illustrates how the real‑imaginary coupling
strength varies with position to steer the field toward the steady
state.  Output: fig2_coupling_matrix.pdf.
"""
import numpy as np
import matplotlib.pyplot as plt

x = np.linspace(2, 100, 500)
K = np.sqrt(2 * x * np.log(x))

plt.rcParams.update({'font.size': 11, 'font.family': 'serif'})
fig, ax1 = plt.subplots(figsize=(10, 5))

color1 = 'tab:blue'
ax1.set_xlabel('x')
ax1.set_ylabel(r'$K(x) = \sqrt{2x\ln x}$', color=color1)
ax1.plot(x, K, color=color1, linewidth=2, label=r'$K(x)$')
ax1.tick_params(axis='y', labelcolor=color1)
ax1.set_xlim([2, 100])
ax1.grid(True, alpha=0.3)

ax2 = ax1.twinx()
color2 = 'tab:red'
ax2.set_ylabel(r'$1/\ln x$', color=color2)
ax2.plot(x, 1/np.log(x), color=color2, linestyle='--', linewidth=2,
         label=r'$1/\ln x$')
ax2.tick_params(axis='y', labelcolor=color2)

# Combine legends from both axes
lines1, labels1 = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
ax1.legend(lines1 + lines2, labels1 + labels2, loc='upper right',
           frameon=False)

plt.title('Coupling matrix K(x) and prime density')
plt.tight_layout()
plt.savefig('fig2_coupling_matrix.pdf', dpi=300, bbox_inches='tight')
plt.show()