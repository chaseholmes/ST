# ============================================================================
# λ₃₀(β) CRITICAL BINDING THRESHOLD
# ============================================================================
#
# Question: Does the m=3 eigenvalue cross zero at a geometric β?
#
# If λ₃₀(β*) = 0 at β* = const × π⁻⁶, then β is not a free parameter
# but a geometric critical coupling.
#
# ============================================================================

import numpy as np
from scipy.linalg import eigh_tridiagonal
from scipy.optimize import brentq
import matplotlib.pyplot as plt

def mu(theta):
    return np.cos(theta)**2 * np.sin(theta)

def K_curvature(theta):
    return 7 - 2 * np.tan(theta)**2

def V_effective(theta, m, delta, beta):
    mu_val = mu(theta)
    mu_sq_reg = mu_val**2 + delta**2
    f_theta = mu_sq_reg / (mu_sq_reg + beta)**2
    centrifugal = m**2 * f_theta
    return centrifugal + K_curvature(theta) / 2

def get_lambda_30(beta, N=8000, eps=0.001):
    """Get the lowest m=3 eigenvalue for given β."""
    delta = 1 / (np.pi**6 * np.sqrt(2))
    theta_max = np.pi/2 - eps
    theta_grid = np.linspace(0, theta_max, N)
    dtheta = theta_grid[1] - theta_grid[0]
    
    V = V_effective(theta_grid[1:-1], m=3, delta=delta, beta=beta)
    diagonal = 2.0 / dtheta**2 + V
    off_diagonal = -np.ones(N - 3) / dtheta**2
    
    eigenvalues, _ = eigh_tridiagonal(
        diagonal, off_diagonal, select='i', select_range=(0, 0)
    )
    return eigenvalues[0]

# ============================================================================
# STEP 1: COARSE SCAN — Find where λ₃₀ crosses zero
# ============================================================================

print("=" * 70)
print("  λ₃₀(β) CRITICAL BINDING THRESHOLD")
print("=" * 70)

delta = 1 / (np.pi**6 * np.sqrt(2))

print(f"\n  δ = 1/(π⁶√2) = {delta:.6e}")
print(f"\n  Coarse scan of λ₃₀(β):")
print(f"  {'β':>14s}  {'λ₃₀':>14s}  {'Bound?':>8s}")
print("  " + "-" * 42)

betas_coarse = np.geomspace(1e-6, 1e-2, 50)
lambdas_coarse = []

for b in betas_coarse:
    lam = get_lambda_30(b)
    lambdas_coarse.append(lam)
    if abs(lam) < 50 or (len(lambdas_coarse) >= 2 and 
        lambdas_coarse[-1] * lambdas_coarse[-2] < 0):
        bound = "YES" if lam < 0 else "no"
        print(f"  {b:>14.6e}  {lam:>14.4f}  {bound:>8s}")

lambdas_coarse = np.array(lambdas_coarse)

# Find zero-crossing brackets
crossings = []
for i in range(len(lambdas_coarse) - 1):
    if lambdas_coarse[i] * lambdas_coarse[i+1] < 0:
        crossings.append((betas_coarse[i], betas_coarse[i+1]))
        print(f"\n  ★ Zero crossing between β = {betas_coarse[i]:.6e} and {betas_coarse[i+1]:.6e}")

# ============================================================================
# STEP 2: PRECISE ZERO-CROSSING via Brent's method
# ============================================================================

print(f"\n" + "=" * 70)
print("  STEP 2: PRECISE ZERO-CROSSING")
print("=" * 70)

for bracket in crossings:
    beta_crit = brentq(get_lambda_30, bracket[0], bracket[1], xtol=1e-14)
    lam_check = get_lambda_30(beta_crit)
    
    print(f"\n  β_crit = {beta_crit:.15e}")
    print(f"  λ₃₀(β_crit) = {lam_check:.6e}  (should be ~0)")

# ============================================================================
# STEP 3: TEST AGAINST GEOMETRIC CANDIDATES
# ============================================================================

print(f"\n" + "=" * 70)
print("  STEP 3: GEOMETRIC IDENTIFICATION OF β_crit")
print("=" * 70)

# Build candidates from the framework's geometric constants
pi = np.pi
sqrt2 = np.sqrt(2)
sqrt3 = np.sqrt(3)
theta_e = np.arccos(np.sqrt(2/3))
mu_peak = 2*sqrt3/9
phi_g = np.sqrt(2/3)  # cos(θ_e)
psi_g = np.sqrt(1/3)  # sin(θ_e)
alpha = 1/137.035999178

candidates = {
    # Powers of π and δ
    'δ = 1/(π⁶√2)':            delta,
    'δ²':                       delta**2,
    'δ² × 2':                   delta**2 * 2,
    'δ² × 3':                   delta**2 * 3,
    'δ² × π':                   delta**2 * pi,
    'δ² × 2π':                  delta**2 * 2*pi,
    '1/(π⁶)':                   1/pi**6,
    '1/(2π⁶)':                  1/(2*pi**6),
    '1/(3π⁶)':                  1/(3*pi**6),
    '1/(4π⁶)':                  1/(4*pi**6),
    '1/(6π⁶)':                  1/(6*pi**6),
    '1/(9π⁶)':                  1/(9*pi**6),
    '1/(12π⁶)':                 1/(12*pi**6),
    '1/(18π⁶)':                 1/(18*pi**6),
    '1/(36π⁶)':                 1/(36*pi**6),
    
    # With √2
    '1/(12π⁶√2)':              1/(12*pi**6*sqrt2),
    '1/(6π⁶√2)':               1/(6*pi**6*sqrt2),
    '1/(9π⁶√2)':               1/(9*pi**6*sqrt2),
    '1/(18π⁶√2)':              1/(18*pi**6*sqrt2),
    
    # From framework constants
    'δ/9':                      delta/9,
    'δ/(3π)':                   delta/(3*pi),
    'δ/(4π)':                   delta/(4*pi),
    'δ/12':                     delta/12,
    'δ/(12π)':                  delta/(12*pi),
    'μ_peak × δ':               mu_peak * delta,
    'μ_peak² × δ':              mu_peak**2 * delta,
    'μ_peak/π⁶':                mu_peak / pi**6,
    
    # Compound expressions
    'δ × ψ_g':                  delta * psi_g,
    'δ × φ_g':                  delta * phi_g,
    'δ × ψ_g²':                delta * psi_g**2,
    'δ × μ_peak':               delta * mu_peak,
    '(∫μdθ)×δ':                (1/3) * delta,
    '(∫μdθ)²×δ':               (1/9) * delta,
    
    # α-based
    'α²':                       alpha**2,
    'α² × π':                   alpha**2 * pi,
    'α² × 2':                   alpha**2 * 2,
    'α/(4π)²':                  alpha/(4*pi)**2,
}

# Sort by closeness to β_crit
scored = []
for name, val in candidates.items():
    if val > 0:
        ratio = beta_crit / val
        log_ratio = abs(np.log(ratio))
        scored.append((log_ratio, ratio, name, val))

scored.sort()

print(f"\n  β_crit = {beta_crit:.10e}")
print(f"\n  {'Candidate':>25s}  {'Value':>14s}  {'β_crit/cand':>14s}  {'Match'}")
print("  " + "-" * 70)

for log_r, ratio, name, val in scored[:25]:
    if 0.1 < ratio < 10:
        marker = " ★★★" if abs(ratio-1) < 0.005 else " ★★" if abs(ratio-1) < 0.02 else " ★" if abs(ratio-1) < 0.05 else ""
        print(f"  {name:>25s}  {val:>14.6e}  {ratio:>14.8f}{marker}")
        
        # Check small integer and simple fraction multiples
        for num in [1, 2, 3, 4, 6, 8, 9, 12]:
            for den in [1, 2, 3, 4, 6, 8, 9, 12]:
                test_ratio = ratio * den / num
                if abs(test_ratio - 1) < 0.003 and not (num == 1 and den == 1):
                    print(f"  {'→ ×'+str(num)+'/'+str(den):>25s}  {val*num/den:>14.6e}  {beta_crit/(val*num/den):>14.8f} ★★★")

# ============================================================================
# STEP 4: RESOLUTION DEPENDENCE — Does β_crit sharpen with N?
# ============================================================================

print(f"\n" + "=" * 70)
print("  STEP 4: RESOLUTION DEPENDENCE OF β_crit")
print("=" * 70)

print(f"\n  Does β_crit converge as N → ∞?")
print(f"\n  {'N':>8s}  {'β_crit':>20s}  {'Δβ/β (vs prev)':>18s}")
print("  " + "-" * 52)

prev_beta = None
N_values = [1000, 2000, 4000, 8000, 16000]

for N in N_values:
    def lam30_N(beta):
        return get_lambda_30(beta, N=N)
    
    # Find bracket first
    b_lo, b_hi = 5e-5, 2e-4
    if lam30_N(b_lo) * lam30_N(b_hi) > 0:
        # Try wider bracket
        b_lo, b_hi = 1e-5, 5e-4
    
    try:
        bc = brentq(lam30_N, b_lo, b_hi, xtol=1e-14)
        if prev_beta is not None:
            rel_change = abs(bc - prev_beta) / bc * 100
            print(f"  {N:>8d}  {bc:>20.12e}  {rel_change:>18.8f}%")
        else:
            print(f"  {N:>8d}  {bc:>20.12e}  {'—':>18s}")
        prev_beta = bc
    except:
        print(f"  {N:>8d}  {'FAILED':>20s}")

beta_converged = prev_beta
print(f"\n  Converged β_crit = {beta_converged:.12e}")

# ============================================================================
# STEP 5: ε DEPENDENCE — Does β_crit depend on the θ cutoff?
# ============================================================================

print(f"\n" + "=" * 70)
print("  STEP 5: CUTOFF DEPENDENCE OF β_crit")
print("=" * 70)

print(f"\n  Does β_crit change with ε (θ cutoff)?")
print(f"\n  {'ε':>12s}  {'β_crit':>20s}  {'β_crit/[1/(12π⁶)]':>20s}")
print("  " + "-" * 58)

ref_val = 1/(12*pi**6)

for eps in [0.01, 0.005, 0.002, 0.001, 0.0005]:
    def lam30_eps(beta):
        return get_lambda_30(beta, N=8000, eps=eps)
    
    try:
        b_lo, b_hi = 1e-5, 1e-2
        # Make sure bracket is valid
        f_lo = lam30_eps(b_lo)
        f_hi = lam30_eps(b_hi)
        if f_lo * f_hi > 0:
            # narrow bracket
            for b_test in np.geomspace(1e-5, 1e-2, 100):
                f_test = lam30_eps(b_test)
                if f_lo * f_test < 0:
                    b_hi = b_test
                    break
                b_lo = b_test
                f_lo = f_test
        
        bc = brentq(lam30_eps, b_lo, b_hi, xtol=1e-14)
        print(f"  {eps:>12.4f}  {bc:>20.12e}  {bc/ref_val:>20.8f}")
    except Exception as e:
        print(f"  {eps:>12.4f}  {'FAILED':>20s}  {str(e)[:30]}")

# ============================================================================
# STEP 6: CHARACTERIZE THE THRESHOLD
# ============================================================================

print(f"\n" + "=" * 70)
print("  STEP 6: PHYSICS AT THE THRESHOLD")
print("=" * 70)

# Compute several quantities at β_crit
delta_val = 1 / (pi**6 * sqrt2)

# Fine scan around β_crit
betas_fine = np.linspace(beta_crit * 0.8, beta_crit * 1.2, 100)
lams_fine = [get_lambda_30(b) for b in betas_fine]

# Slope dλ/dβ at threshold
db = beta_crit * 1e-6
slope = (get_lambda_30(beta_crit + db) - get_lambda_30(beta_crit - db)) / (2*db)

print(f"""
  AT THE CRITICAL THRESHOLD:
  ──────────────────────────
  β_crit = {beta_crit:.10e}
  λ₃₀(β_crit) ≈ 0
  
  dλ₃₀/dβ at threshold = {slope:.4f}
  
  GEOMETRIC COMPARISONS:
  ──────────────────────
  δ = 1/(π⁶√2)     = {delta_val:.10e}
  β_crit/δ          = {beta_crit/delta_val:.10f}
  
  1/(12π⁶)          = {ref_val:.10e}
  β_crit/[1/(12π⁶)] = {beta_crit/ref_val:.10f}
  
  δ/9               = {delta_val/9:.10e}
  β_crit/(δ/9)      = {beta_crit/(delta_val/9):.10f}
  
  δ²×2π⁶            = {delta_val**2 * 2 * pi**6:.10e}
  
  Key ratios:
  β_crit × π⁶       = {beta_crit * pi**6:.10f}
  β_crit × 12π⁶     = {beta_crit * 12 * pi**6:.10f}
  β_crit × π⁶√2     = {beta_crit * pi**6 * sqrt2:.10f}
  β_crit / δ         = {beta_crit / delta_val:.10f}
  β_crit × 9/δ       = {beta_crit * 9 / delta_val:.10f}
""")

# ============================================================================
# VISUALIZATION
# ============================================================================

fig, axes = plt.subplots(1, 3, figsize=(18, 5))
fig.suptitle('λ₃₀(β): Critical Binding Threshold', fontsize=14, fontweight='bold')

# Panel 1: Full range (log scale)
ax = axes[0]
mask = np.array(lambdas_coarse) > -2000  # avoid extreme values
ax.semilogx(betas_coarse[mask], lambdas_coarse[mask], 'b.-', linewidth=2, markersize=4)
ax.axhline(y=0, color='red', linewidth=2, linestyle='--', label='λ₃₀ = 0 (threshold)')
ax.axvline(x=beta_crit, color='purple', linewidth=2, linestyle='--', label=f'β_crit = {beta_crit:.2e}')
ax.axvline(x=ref_val, color='green', linewidth=1.5, linestyle=':', label=f'1/(12π⁶) = {ref_val:.2e}')
ax.set_xlabel('β')
ax.set_ylabel('λ₃₀')
ax.set_title('Full Range')
ax.legend(fontsize=8)
ax.grid(True, alpha=0.3)
ax.set_ylim(-500, 500)

# Panel 2: Zoom near threshold
ax = axes[1]
ax.plot(betas_fine, lams_fine, 'b-', linewidth=2)
ax.axhline(y=0, color='red', linewidth=2, linestyle='--')
ax.axvline(x=beta_crit, color='purple', linewidth=2, linestyle='--', label=f'β_crit')
ax.axvline(x=ref_val, color='green', linewidth=1.5, linestyle=':', label=f'1/(12π⁶)')
ax.fill_between(betas_fine, lams_fine, 0, where=np.array(lams_fine)<0, alpha=0.2, color='blue', label='Bound (λ<0)')
ax.fill_between(betas_fine, lams_fine, 0, where=np.array(lams_fine)>0, alpha=0.2, color='red', label='Unbound (λ>0)')
ax.set_xlabel('β')
ax.set_ylabel('λ₃₀')
ax.set_title('Near Threshold')
ax.legend(fontsize=8)
ax.grid(True, alpha=0.3)

# Panel 3: Convergence with N
ax = axes[2]
N_conv = [1000, 2000, 4000, 8000, 16000]
bc_conv = []
for N in N_conv:
    def lam_N(beta):
        return get_lambda_30(beta, N=N)
    try:
        bc = brentq(lam_N, 5e-5, 2e-4, xtol=1e-14)
        bc_conv.append(bc)
    except:
        bc_conv.append(np.nan)

ax.semilogx(N_conv, bc_conv, 'ko-', linewidth=2, markersize=8)
ax.axhline(y=ref_val, color='green', linewidth=2, linestyle='--', label=f'1/(12π⁶) = {ref_val:.4e}')
ax.set_xlabel('Grid points N')
ax.set_ylabel('β_crit')
ax.set_title('Convergence of β_crit with N')
ax.legend(fontsize=9)
ax.grid(True, alpha=0.3)

# Add percentage annotation
if len(bc_conv) >= 2 and not np.isnan(bc_conv[-1]):
    pct = (bc_conv[-1] - ref_val) / ref_val * 100
    ax.text(0.5, 0.15, f'β_crit/[1/(12π⁶)] = {bc_conv[-1]/ref_val:.4f}\n({pct:+.2f}%)',
            transform=ax.transAxes, fontsize=11, fontweight='bold',
            ha='center', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.7))

plt.tight_layout()
plt.savefig('/home/claude/beta_threshold.png', dpi=150, bbox_inches='tight')
print(f"  [Plot saved to 'beta_threshold.png']")

# ============================================================================
# VERDICT
# ============================================================================

print("\n" + "=" * 70)
print("  VERDICT")
print("=" * 70)

ratio_12pi6 = beta_converged / ref_val
print(f"""
  β_crit (converged) = {beta_converged:.10e}
  1/(12π⁶)           = {ref_val:.10e}
  
  Ratio: {ratio_12pi6:.8f}  ({(ratio_12pi6-1)*100:+.4f}%)
  
  β_crit × π⁶ = {beta_converged * pi**6:.8f}  (if exactly 1/12, this = 0.08333...)
  1/12 = {1/12:.8f}
""")

if abs(ratio_12pi6 - 1) < 0.1:
    print(f"  β_crit ≈ 1/(12π⁶) within {abs(ratio_12pi6-1)*100:.1f}%")
    print(f"  This is suggestive but not exact.")
    
    # Check nearby simple fractions
    for num, den in [(1,12), (1,11), (1,13), (1,9), (1,10), (1,8),
                     (1,6), (2,9), (2,11), (3,8)]:
        test = num/(den * pi**6)
        r = beta_converged / test
        if abs(r - 1) < 0.03:
            print(f"  ★ β_crit ≈ {num}/{den}π⁶ = {test:.6e}, ratio = {r:.6f} ({(r-1)*100:+.3f}%)")

print("\n" + "=" * 70)
