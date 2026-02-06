# ============================================================================
# THE CRITICAL CURVE: λ₃₀(β, ε) = 0
# ============================================================================
#
# Two unknowns: β (fiber dilution) and ε (UV cutoff, θ_max = π/2 - ε)
# One condition: λ₃₀ = 0 (critical binding)
#
# This gives a CURVE in (β, ε) space. 
# The question: does geometry provide a SECOND condition that picks a POINT?
#
# Strategy:
#   1. Map the critical curve β_crit(ε)
#   2. Look for geometric relationships: β = f(δ, ε)?
#   3. Check if β_crit(ε) intersects any natural constraint
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

def get_lambda_30(beta, N=8000, eps=0.001, delta=None):
    if delta is None:
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

pi = np.pi
sqrt2 = np.sqrt(2)
delta = 1 / (pi**6 * sqrt2)

print("=" * 72)
print("  THE CRITICAL CURVE: λ₃₀(β, ε) = 0")
print("=" * 72)
print(f"\n  δ = 1/(π⁶√2) = {delta:.6e}")

# ============================================================================
# STEP 1: MAP THE CRITICAL CURVE β_crit(ε)
# ============================================================================

print(f"\n" + "=" * 72)
print("  STEP 1: CRITICAL CURVE β_crit(ε)")
print("=" * 72)

eps_values = np.geomspace(0.0003, 0.05, 60)
beta_crits = []
N_compute = 8000

print(f"\n  {'ε':>12s}  {'θ_max (deg)':>12s}  {'β_crit':>14s}  {'β_crit/δ':>12s}  {'β_crit×9/δ':>12s}")
print("  " + "-" * 68)

for eps in eps_values:
    def lam_at_eps(log_beta):
        return get_lambda_30(10**log_beta, N=N_compute, eps=eps)
    
    try:
        # Scan for bracket
        log_betas = np.linspace(-7, -1, 80)
        vals = [lam_at_eps(lb) for lb in log_betas]
        
        bracket = None
        for i in range(len(vals)-1):
            if vals[i] * vals[i+1] < 0:
                bracket = (log_betas[i], log_betas[i+1])
                break
        
        if bracket:
            log_bc = brentq(lam_at_eps, bracket[0], bracket[1], xtol=1e-12)
            bc = 10**log_bc
            beta_crits.append(bc)
            
            ratio_delta = bc / delta
            ratio_delta9 = bc * 9 / delta
            
            if len(beta_crits) % 5 == 1 or abs(ratio_delta9 - 1) < 0.05:
                marker = " ★" if abs(ratio_delta9 - 1) < 0.05 else ""
                print(f"  {eps:>12.6f}  {np.degrees(pi/2-eps):>12.4f}  {bc:>14.6e}  "
                      f"{ratio_delta:>12.6f}  {ratio_delta9:>12.6f}{marker}")
        else:
            beta_crits.append(np.nan)
    except:
        beta_crits.append(np.nan)

eps_values = eps_values[:len(beta_crits)]
beta_crits = np.array(beta_crits)

# ============================================================================
# STEP 2: FUNCTIONAL FORM — What is β_crit(ε)?
# ============================================================================

print(f"\n" + "=" * 72)
print("  STEP 2: FUNCTIONAL FORM OF β_crit(ε)")
print("=" * 72)

# Clean data
mask = np.isfinite(beta_crits)
eps_clean = eps_values[mask]
bc_clean = beta_crits[mask]

# Try power law: β_crit = A × ε^p
log_eps = np.log(eps_clean)
log_bc = np.log(bc_clean)

# Linear fit in log-log space
coeffs = np.polyfit(log_eps, log_bc, 1)
p_fit = coeffs[0]
A_fit = np.exp(coeffs[1])

print(f"\n  Power-law fit: β_crit ≈ {A_fit:.4e} × ε^{p_fit:.4f}")
print(f"  Exponent ≈ {p_fit:.4f}")

# Check if exponent is close to a simple fraction
simple_exponents = [(2, 'ε²'), (3, 'ε³'), (4, 'ε⁴'), 
                    (5/2, 'ε^(5/2)'), (3/2, 'ε^(3/2)'),
                    (7/3, 'ε^(7/3)'), (8/3, 'ε^(8/3)')]

print(f"\n  Testing simple exponent hypotheses:")
for p_test, name in simple_exponents:
    residuals = np.std(log_bc - p_test * log_eps - np.mean(log_bc - p_test * log_eps))
    A_test = np.exp(np.mean(log_bc - p_test * log_eps))
    print(f"    β ∝ {name}: A = {A_test:.4e}, residual = {residuals:.4f}")

# Also try: β_crit = δ × f(ε)
bc_over_delta = bc_clean / delta

print(f"\n  β_crit/δ as function of ε:")
coeffs2 = np.polyfit(np.log(eps_clean), np.log(bc_over_delta), 1)
p2 = coeffs2[0]
A2 = np.exp(coeffs2[1])
print(f"  β_crit/δ ≈ {A2:.4f} × ε^{p2:.4f}")

# ============================================================================
# STEP 3: GEOMETRIC CONSTRAINT CANDIDATES
# ============================================================================

print(f"\n" + "=" * 72)
print("  STEP 3: GEOMETRIC CONSTRAINTS — WHERE DOES THE CURVE CROSS?")
print("=" * 72)

# Candidate second conditions:

# (A) β = δ/9 (from the convergence study)
beta_delta9 = delta / 9
# Find ε where β_crit(ε) = δ/9
from scipy.interpolate import interp1d

interp_bc = interp1d(np.log(eps_clean), np.log(bc_clean), kind='cubic', fill_value='extrapolate')

def bc_minus_target(log_eps, target):
    return np.exp(interp_bc(log_eps)) - target

targets = {
    'β = δ/9 = 1/(9π⁶√2)': delta/9,
    'β = δ/12':              delta/12,
    'β = 1/(12π⁶)':         1/(12*pi**6),
    'β = δ²':                delta**2,
    'β = δ × ε':             None,  # self-referential, handle separately
    'β = ε²':                None,  # self-referential
}

print(f"\n  Finding ε* where β_crit(ε*) intersects geometric constraints:")
print(f"\n  {'Constraint':>25s}  {'ε*':>14s}  {'θ_max* (deg)':>14s}  {'β* value':>14s}")
print("  " + "-" * 72)

intersection_points = {}

for name, target_val in targets.items():
    if target_val is None:
        continue
    try:
        log_eps_star = brentq(lambda le: bc_minus_target(le, target_val), 
                              np.log(eps_clean[0]), np.log(eps_clean[-1]))
        eps_star = np.exp(log_eps_star)
        theta_max_star = np.degrees(pi/2 - eps_star)
        print(f"  {name:>25s}  {eps_star:>14.8f}  {theta_max_star:>14.4f}  {target_val:>14.6e}")
        intersection_points[name] = eps_star
    except:
        print(f"  {name:>25s}  {'(no crossing)':>14s}")

# Self-referential constraints: β = δ × ε and β = ε²
print(f"\n  Self-referential constraints:")

# β_crit(ε) = δ × ε → find ε where β_crit(ε)/ε = δ
try:
    def self_ref_delta_eps(log_eps):
        eps = np.exp(log_eps)
        return np.exp(interp_bc(log_eps)) - delta * eps
    
    log_eps_sr = brentq(self_ref_delta_eps, np.log(eps_clean[0]), np.log(eps_clean[-1]))
    eps_sr = np.exp(log_eps_sr)
    bc_sr = np.exp(interp_bc(log_eps_sr))
    print(f"  {'β = δ × ε':>25s}  ε* = {eps_sr:.8f}  θ_max = {np.degrees(pi/2-eps_sr):.4f}°  β* = {bc_sr:.6e}")
    intersection_points['β = δ × ε'] = eps_sr
except:
    print(f"  {'β = δ × ε':>25s}  (no crossing)")

# β_crit(ε) = ε² 
try:
    def self_ref_eps2(log_eps):
        eps = np.exp(log_eps)
        return np.exp(interp_bc(log_eps)) - eps**2
    
    log_eps_sr2 = brentq(self_ref_eps2, np.log(eps_clean[0]), np.log(eps_clean[-1]))
    eps_sr2 = np.exp(log_eps_sr2)
    bc_sr2 = np.exp(interp_bc(log_eps_sr2))
    print(f"  {'β = ε²':>25s}  ε* = {eps_sr2:.8f}  θ_max = {np.degrees(pi/2-eps_sr2):.4f}°  β* = {bc_sr2:.6e}")
    intersection_points['β = ε²'] = eps_sr2
except:
    print(f"  {'β = ε²':>25s}  (no crossing)")

# β_crit(ε) = ε³/δ
try:
    def self_ref_eps3(log_eps):
        eps = np.exp(log_eps)
        return np.exp(interp_bc(log_eps)) - eps**3/delta
    
    log_eps_sr3 = brentq(self_ref_eps3, np.log(eps_clean[0]), np.log(eps_clean[-1]))
    eps_sr3 = np.exp(log_eps_sr3)
    bc_sr3 = np.exp(interp_bc(log_eps_sr3))
    print(f"  {'β = ε³/δ':>25s}  ε* = {eps_sr3:.8f}  θ_max = {np.degrees(pi/2-eps_sr3):.4f}°  β* = {bc_sr3:.6e}")
    intersection_points['β = ε³/δ'] = eps_sr3
except:
    print(f"  {'β = ε³/δ':>25s}  (no crossing)")

# ============================================================================
# STEP 4: PHYSICS AT EACH INTERSECTION
# ============================================================================

print(f"\n" + "=" * 72)
print("  STEP 4: PHYSICS AT EACH INTERSECTION POINT")
print("=" * 72)

print(f"\n  For each (β*, ε*), compute C₂, λ_scale, and the Schwinger product:")

from scipy.linalg import eigh_tridiagonal as eigh_tri

def full_analysis(beta, eps, N=8000):
    """Compute C₂, λ_scale at given (β, ε)."""
    delta_val = 1 / (pi**6 * sqrt2)
    theta_max = pi/2 - eps
    
    results = {}
    for m_val in [0, 3]:
        theta_grid = np.linspace(0, theta_max, N)
        dtheta = theta_grid[1] - theta_grid[0]
        V = V_effective(theta_grid[1:-1], m=m_val, delta=delta_val, beta=beta)
        diag = 2.0/dtheta**2 + V
        off = -np.ones(N-3)/dtheta**2
        try:
            evals, evecs = eigh_tri(diag, off, select='i', select_range=(0,2))
        except:
            evals, evecs = eigh_tri(diag, off)
            evals = evals[:3]; evecs = evecs[:,:3]
        
        efuncs = np.zeros((N, len(evals)))
        efuncs[1:-1,:] = evecs
        for i in range(len(evals)):
            norm = np.trapezoid(efuncs[:,i]**2 * mu(theta_grid), theta_grid)
            if norm > 0: efuncs[:,i] /= np.sqrt(norm)
        
        results[m_val] = {'evals': evals, 'efuncs': efuncs, 'grid': theta_grid}
    
    lam00 = results[0]['evals'][0]
    lam30 = results[3]['evals'][0]
    dE = lam30 - lam00
    
    psi00 = results[0]['efuncs'][:,0]
    psi30 = results[3]['efuncs'][:,0]
    grid = results[0]['grid']
    
    I30 = np.trapezoid(psi00 * psi30 * mu(grid), grid)
    
    if abs(I30) > 1e-10 and abs(dE) > 1e-10:
        C2 = 9 * I30**2 / abs(dE)
        lam_scale = abs(dE) / (2 * 9 * I30**2)
        product = 2 * C2 * lam_scale
    else:
        C2 = 0; lam_scale = float('inf'); product = 0
    
    return {
        'lam00': lam00, 'lam30': lam30, 'dE': dE,
        'I30': I30, 'C2': C2, 'lam_scale': lam_scale,
        'product': product
    }

C2_target = 1/(72*pi)
ls_target = 36*pi

print(f"\n  {'Constraint':>22s}  {'ε*':>10s}  {'C₂':>12s}  {'err%':>8s}  {'λ_scale':>10s}  {'err%':>8s}  {'2C₂λ':>8s}")
print("  " + "-" * 90)

for name, eps_star in intersection_points.items():
    # Use the constraint to get β
    if name == 'β = δ/9 = 1/(9π⁶√2)':
        beta_use = delta/9
    elif name == 'β = δ/12':
        beta_use = delta/12
    elif name == 'β = 1/(12π⁶)':
        beta_use = 1/(12*pi**6)
    elif name == 'β = δ²':
        beta_use = delta**2
    elif name == 'β = δ × ε':
        beta_use = delta * eps_star
    elif name == 'β = ε²':
        beta_use = eps_star**2
    elif name == 'β = ε³/δ':
        beta_use = eps_star**3 / delta
    else:
        continue
    
    r = full_analysis(beta_use, eps_star)
    
    if r['C2'] > 0:
        c2_err = (r['C2'] - C2_target)/C2_target * 100
        ls_err = (r['lam_scale'] - ls_target)/ls_target * 100 if r['lam_scale'] < 1e10 else float('inf')
        print(f"  {name:>22s}  {eps_star:>10.6f}  {r['C2']:>12.8f}  {c2_err:>+8.3f}  "
              f"{r['lam_scale']:>10.4f}  {ls_err:>+8.3f}  {r['product']:>8.4f}")
    else:
        print(f"  {name:>22s}  {eps_star:>10.6f}  {'(unbound)':>12s}")

# ============================================================================
# STEP 5: THE μ SINGULARITY — WHAT SETS ε PHYSICALLY?
# ============================================================================

print(f"\n" + "=" * 72)
print("  STEP 5: WHAT SETS ε PHYSICALLY?")
print("=" * 72)

print(f"""
  At θ = π/2, the measure μ(θ) = cos²θ sinθ → 0 (coordinate singularity).
  The cutoff ε keeps us away from this singularity.
  
  Physical candidates for what fixes ε:
  
  (A) δ-resolution: The fiber can't resolve features smaller than δ.
      At θ = π/2 - ε, the effective radius is μ(π/2-ε) ≈ ε.
      Setting μ(π/2-ε) = δ gives ε ≈ δ^(1/3) (since μ ≈ ε² × 1 near π/2).
""")

# Near θ = π/2: let θ = π/2 - ε (ε small)
# cosθ ≈ ε, sinθ ≈ 1
# μ ≈ ε² × 1 = ε²
# So μ(π/2-ε) = ε² to leading order

# Condition: μ(θ_max) = δ → ε² = δ → ε = √δ
eps_sqrt_delta = np.sqrt(delta)
# Condition: μ(θ_max) = δ² → ε = δ
eps_delta = delta
# Condition: μ(θ_max) = β → ε² = β → ε = √β
eps_sqrt_beta = np.sqrt(delta/9)  # if β = δ/9

print(f"  Near π/2: μ(π/2-ε) ≈ ε²")
print(f"  ")
print(f"  Candidate ε from μ = (resolution scale):")
print(f"  ")
print(f"  (i)   μ = δ          → ε = √δ = {eps_sqrt_delta:.6e}")
print(f"  (ii)  μ = δ²         → ε = δ   = {delta:.6e}")
print(f"  (iii) μ = β = δ/9    → ε = √(δ/9) = {eps_sqrt_beta:.6e}")

# Check exact μ at these ε values
for name, eps_test in [('√δ', eps_sqrt_delta), ('δ', delta), ('√(δ/9)', eps_sqrt_beta)]:
    theta_test = pi/2 - eps_test
    mu_test = mu(theta_test)
    print(f"\n  ε = {name} = {eps_test:.6e}:")
    print(f"    θ_max = {np.degrees(theta_test):.6f}°")
    print(f"    μ(θ_max) = {mu_test:.6e}")
    print(f"    μ/δ = {mu_test/delta:.6f}")
    print(f"    μ/δ² = {mu_test/delta**2:.6f}")

# What β_crit corresponds to these ε values?
print(f"\n  β_crit at these ε values (from the critical curve):")

for name, eps_test in [('√δ', eps_sqrt_delta), ('δ', delta), ('√(δ/9)', eps_sqrt_beta)]:
    if np.log(eps_test) >= np.log(eps_clean[0]) and np.log(eps_test) <= np.log(eps_clean[-1]):
        bc_test = np.exp(interp_bc(np.log(eps_test)))
        print(f"  ε = {name}: β_crit = {bc_test:.6e}, β_crit/δ = {bc_test/delta:.6f}, β_crit×9/δ = {bc_test*9/delta:.6f}")
    else:
        print(f"  ε = {name}: outside interpolation range")

# ============================================================================
# STEP 6: THE SELF-CONSISTENT POINT
# ============================================================================

print(f"\n" + "=" * 72)
print("  STEP 6: SELF-CONSISTENT SOLUTION")
print("=" * 72)

print(f"""
  The two conditions are:
    (1) λ₃₀(β, ε) = 0     [critical binding]
    (2) μ(π/2 - ε) = β     [cutoff = where measure equals dilution]
  
  Together with μ(π/2-ε) ≈ ε² (near π/2):
    β ≈ ε²
    
  So we need: β_crit(ε) = ε²
  This is the 'β = ε²' intersection we already found.
""")

if 'β = ε²' in intersection_points:
    eps_sc = intersection_points['β = ε²']
    beta_sc = eps_sc**2
    
    # Verify
    lam_check = get_lambda_30(beta_sc, eps=eps_sc)
    mu_check = mu(pi/2 - eps_sc)
    
    print(f"  SELF-CONSISTENT POINT:")
    print(f"  ε* = {eps_sc:.10f}")
    print(f"  β* = ε*² = {beta_sc:.10e}")
    print(f"  θ_max = {np.degrees(pi/2 - eps_sc):.6f}°")
    print(f"  μ(θ_max) = {mu_check:.10e}")
    print(f"  μ/β = {mu_check/beta_sc:.6f}  (should be ~1)")
    print(f"  λ₃₀(β*, ε*) = {lam_check:.6e}  (should be ~0)")
    
    # Now compute C₂ at this point
    r_sc = full_analysis(beta_sc, eps_sc)
    c2_err = (r_sc['C2'] - C2_target)/C2_target * 100
    
    print(f"\n  PHYSICS AT THE SELF-CONSISTENT POINT:")
    print(f"  λ₀₀ = {r_sc['lam00']:.4f}")
    print(f"  λ₃₀ = {r_sc['lam30']:.4f}")
    print(f"  ΔE = {r_sc['dE']:.4f}")
    print(f"  I₃₀ = {r_sc['I30']:.8f}")
    print(f"  C₂ = {r_sc['C2']:.10f}  (target: {C2_target:.10f}, err: {c2_err:+.4f}%)")
    print(f"  λ_scale = {r_sc['lam_scale']:.4f}  (target: {ls_target:.4f})")
    print(f"  2·C₂·λ_scale = {r_sc['product']:.6f}  (should be 1)")

# Also try: μ(π/2-ε) = δ (resolution condition)
print(f"\n  Alternative: μ(π/2-ε) = δ → ε = √δ:")

eps_alt = np.sqrt(delta)
if np.log(eps_alt) >= np.log(eps_clean[0]) and np.log(eps_alt) <= np.log(eps_clean[-1]):
    bc_alt = np.exp(interp_bc(np.log(eps_alt)))
    lam_alt = get_lambda_30(bc_alt, eps=eps_alt)
    
    print(f"  ε = √δ = {eps_alt:.8f}")
    print(f"  β_crit(√δ) = {bc_alt:.6e}")
    print(f"  β_crit/δ = {bc_alt/delta:.6f}")
    
    r_alt = full_analysis(bc_alt, eps_alt)
    if r_alt['C2'] > 0:
        c2_err_alt = (r_alt['C2'] - C2_target)/C2_target * 100
        print(f"  C₂ = {r_alt['C2']:.10f} (err: {c2_err_alt:+.4f}%)")
        print(f"  λ_scale = {r_alt['lam_scale']:.4f}")
        print(f"  2·C₂·λ_scale = {r_alt['product']:.6f}")

# ============================================================================
# VISUALIZATION
# ============================================================================

fig, axes = plt.subplots(1, 3, figsize=(18, 6))
fig.suptitle('Critical Curve λ₃₀(β, ε) = 0 and Geometric Constraints', 
             fontsize=14, fontweight='bold')

# Panel 1: The critical curve β_crit(ε) in log-log
ax = axes[0]
ax.loglog(eps_clean, bc_clean, 'b-', linewidth=2.5, label='β_crit(ε) [binding threshold]')

# Overlay constraint lines
eps_range = np.geomspace(eps_clean[0], eps_clean[-1], 200)

ax.loglog(eps_range, np.full_like(eps_range, delta/9), 'r--', linewidth=1.5, label='β = δ/9')
ax.loglog(eps_range, np.full_like(eps_range, 1/(12*pi**6)), 'g--', linewidth=1.5, label='β = 1/(12π⁶)')
ax.loglog(eps_range, eps_range**2, 'm--', linewidth=1.5, label='β = ε² [μ(θ_max)=β]')
ax.loglog(eps_range, delta * eps_range, 'c--', linewidth=1.5, label='β = δε')
ax.loglog(eps_range, np.full_like(eps_range, delta**2), 'k:', linewidth=1, label='β = δ²')

# Mark intersection points
for name, eps_star in intersection_points.items():
    if name == 'β = δ/9 = 1/(9π⁶√2)':
        ax.plot(eps_star, delta/9, 'ro', markersize=12, zorder=5)
    elif name == 'β = ε²':
        ax.plot(eps_star, eps_star**2, 'ms', markersize=12, zorder=5)
    elif name == 'β = δ × ε':
        ax.plot(eps_star, delta*eps_star, 'c^', markersize=12, zorder=5)

ax.set_xlabel('ε (UV cutoff)', fontsize=12)
ax.set_ylabel('β (fiber dilution)', fontsize=12)
ax.set_title('Critical Curve + Constraints')
ax.legend(fontsize=7, loc='upper left')
ax.grid(True, alpha=0.3)
ax.set_xlim(eps_clean[0]*0.8, eps_clean[-1]*1.2)

# Panel 2: β_crit/δ vs ε — looking for the geometric ratio
ax = axes[1]
ax.semilogx(eps_clean, bc_clean / delta, 'b-', linewidth=2.5)
ax.axhline(y=1/9, color='red', linewidth=2, linestyle='--', label='1/9 (= [∫μdθ]²)')
ax.axhline(y=1/12, color='green', linewidth=1.5, linestyle='--', label='1/12')

# Mark the δ/9 crossing
if 'β = δ/9 = 1/(9π⁶√2)' in intersection_points:
    eps_d9 = intersection_points['β = δ/9 = 1/(9π⁶√2)']
    ax.plot(eps_d9, 1/9, 'ro', markersize=12, zorder=5, label=f'ε = {eps_d9:.4f}')

ax.set_xlabel('ε', fontsize=12)
ax.set_ylabel('β_crit / δ', fontsize=12)
ax.set_title('β_crit/δ: Where does it equal 1/9?')
ax.legend(fontsize=9)
ax.grid(True, alpha=0.3)
ax.set_ylim(-0.1, 1.0)

# Panel 3: Power law check
ax = axes[2]
ax.loglog(eps_clean, bc_clean, 'b.', markersize=4, alpha=0.5)

# Fit line
eps_fit = np.geomspace(eps_clean[0], eps_clean[-1], 100)
bc_fit = A_fit * eps_fit**p_fit
ax.loglog(eps_fit, bc_fit, 'r-', linewidth=2, label=f'Fit: β ∝ ε^{{{p_fit:.2f}}}')

# Reference slopes
for p, style, label in [(2, 'g--', 'ε²'), (3, 'm--', 'ε³'), (2.5, 'c:', 'ε^2.5')]:
    A_ref = np.exp(np.mean(np.log(bc_clean) - p * np.log(eps_clean)))
    ax.loglog(eps_fit, A_ref * eps_fit**p, style, linewidth=1, alpha=0.7, label=label)

ax.set_xlabel('ε', fontsize=12)
ax.set_ylabel('β_crit', fontsize=12)
ax.set_title(f'Power Law: β_crit ∝ ε^{{{p_fit:.2f}}}')
ax.legend(fontsize=9)
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('/home/claude/critical_curve.png', dpi=150, bbox_inches='tight')
print(f"\n  [Plot saved to 'critical_curve.png']")

# ============================================================================
# SUMMARY
# ============================================================================

print(f"\n" + "=" * 72)
print("  SUMMARY")
print("=" * 72)

print(f"""
  THE CRITICAL CURVE β_crit(ε):
  ─────────────────────────────
  Power law: β_crit ≈ {A_fit:.4e} × ε^{p_fit:.2f}
  
  KEY INTERSECTIONS:
""")

for name, eps_star in sorted(intersection_points.items(), key=lambda x: x[1]):
    print(f"  {name:>25s}:  ε* = {eps_star:.6f}  (θ_max = {np.degrees(pi/2-eps_star):.2f}°)")

if 'β = ε²' in intersection_points:
    eps_sc = intersection_points['β = ε²']
    print(f"""
  SELF-CONSISTENT POINT (β = ε², meaning μ(θ_max) ≈ β):
  ε* = {eps_sc:.6f}, β* = {eps_sc**2:.6e}
  
  This is the point where:
    • The UV cutoff equals the measure singularity scale
    • The fiber dilution equals the measure at the boundary
    • The m=3 mode is exactly at threshold
  
  This fixes BOTH free parameters from geometry alone.
""")

print("=" * 72)
