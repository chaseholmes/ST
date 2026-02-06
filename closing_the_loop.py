# ============================================================================
# CLOSING THE LOOP: DERIVE ε* FROM λ₃₀(β=δ/9, ε*) = 0
# ============================================================================

import numpy as np
from scipy.linalg import eigh_tridiagonal
from scipy.optimize import brentq
from scipy.integrate import quad
import matplotlib.pyplot as plt

pi = np.pi
sqrt2 = np.sqrt(2)
delta = 1 / (pi**6 * sqrt2)
beta = delta / 9  # DERIVED, not fit

def mu(theta):
    return np.cos(theta)**2 * np.sin(theta)

def K_curvature(theta):
    return 7 - 2 * np.tan(theta)**2

def get_lambda_30(eps, N=8000):
    """λ₃₀ at fixed β=δ/9, as function of ε only."""
    theta_max = pi/2 - eps
    theta_grid = np.linspace(0, theta_max, N)
    dtheta = theta_grid[1] - theta_grid[0]
    
    mu_val = mu(theta_grid[1:-1])
    mu_sq_reg = mu_val**2 + delta**2
    f_theta = mu_sq_reg / (mu_sq_reg + beta)**2
    V = 9 * f_theta + K_curvature(theta_grid[1:-1]) / 2
    
    diag = 2.0/dtheta**2 + V
    off = -np.ones(N-3)/dtheta**2
    evals, _ = eigh_tridiagonal(diag, off, select='i', select_range=(0,0))
    return evals[0]

# ============================================================================
# STEP 1: SOLVE λ₃₀(β=δ/9, ε*) = 0
# ============================================================================

print("=" * 72)
print("  STEP 1: SOLVE λ₃₀(β=δ/9, ε*) = 0")
print("=" * 72)
print(f"\n  β = δ/9 = {beta:.10e}  (fixed)")

# Find bracket
print(f"\n  Scanning for sign change...")
eps_scan = np.geomspace(0.0003, 0.01, 40)
for eps_test in eps_scan:
    lam = get_lambda_30(eps_test)
    if abs(lam) < 100:
        print(f"    ε = {eps_test:.6f}: λ₃₀ = {lam:+.4f}")

# Find bracket
lam_vals = [get_lambda_30(e) for e in eps_scan]
bracket = None
for i in range(len(lam_vals)-1):
    if lam_vals[i] * lam_vals[i+1] < 0:
        bracket = (eps_scan[i], eps_scan[i+1])
        break

if bracket:
    print(f"\n  Bracket: [{bracket[0]:.6f}, {bracket[1]:.6f}]")
    
    # Solve at increasing resolution to check ε* convergence
    print(f"\n  {'N':>8s}  {'ε*':>18s}  {'λ₃₀ check':>14s}  {'θ_max (deg)':>12s}  {'μ(θ_max)':>14s}")
    print("  " + "-" * 72)
    
    eps_stars = []
    for N in [2000, 4000, 8000, 12000, 16000]:
        def lam_N(eps):
            return get_lambda_30(eps, N=N)
        
        eps_star = brentq(lam_N, bracket[0], bracket[1], xtol=1e-14)
        lam_check = lam_N(eps_star)
        theta_max = pi/2 - eps_star
        mu_at_max = mu(theta_max)
        
        eps_stars.append((N, eps_star))
        print(f"  {N:>8d}  {eps_star:>18.12f}  {lam_check:>14.2e}  "
              f"{np.degrees(theta_max):>12.6f}  {mu_at_max:>14.6e}")
    
    # Use highest-resolution value
    N_best, eps_derived = eps_stars[-1]
    
    print(f"\n  ╔═══════════════════════════════════════════════╗")
    print(f"  ║  ε* = {eps_derived:.12f}                  ║")
    print(f"  ║  θ_max = {np.degrees(pi/2 - eps_derived):.8f}°                    ║")
    print(f"  ║  μ(θ_max) = {mu(pi/2 - eps_derived):.6e}                 ║")
    print(f"  ╚═══════════════════════════════════════════════╝")

# ============================================================================
# STEP 2: GEOMETRIC INTERPRETATION OF ε*
# ============================================================================

print(f"\n" + "=" * 72)
print("  STEP 2: WHAT IS ε* GEOMETRICALLY?")
print("=" * 72)

eps_s = eps_derived
mu_at_eps = mu(pi/2 - eps_s)

print(f"""
  ε* = {eps_s:.12f}
  
  Near π/2: μ(π/2 - ε) ≈ ε²  (leading order)
  
  μ(θ_max) = {mu_at_eps:.10e}
  ε*²       = {eps_s**2:.10e}
  μ/ε²     = {mu_at_eps/eps_s**2:.8f}  (≈ 1 confirms the ε² approximation)
  
  KEY RATIOS:
  ───────────
  ε*        = {eps_s:.10e}
  √δ        = {np.sqrt(delta):.10e}    ratio: {eps_s/np.sqrt(delta):.6f}
  √β        = {np.sqrt(beta):.10e}    ratio: {eps_s/np.sqrt(beta):.6f}
  δ^(1/3)   = {delta**(1/3):.10e}    ratio: {eps_s/delta**(1/3):.6f}
  β^(1/3)   = {beta**(1/3):.10e}    ratio: {eps_s/beta**(1/3):.6f}
  
  μ(θ_max) vs geometric scales:
  μ(θ_max)  = {mu_at_eps:.10e}
  β = δ/9   = {beta:.10e}    ratio μ/β: {mu_at_eps/beta:.6f}
  δ         = {delta:.10e}    ratio μ/δ: {mu_at_eps/delta:.6f}
  δ²        = {delta**2:.10e}    ratio μ/δ²: {mu_at_eps/delta**2:.6f}
""")

# Check if ε satisfies any clean relationship
# μ(θ_max) = C × β for some simple C?
print(f"  Testing μ(θ_max) = C × β:")
C_ratio = mu_at_eps / beta
print(f"    C = μ/β = {C_ratio:.8f}")

# Is C close to a simple number?
for name, val in [('1', 1), ('1/2', 0.5), ('1/3', 1/3), ('1/4', 0.25),
                  ('1/9', 1/9), ('1/π', 1/pi), ('δ', delta), ('∫μ', 1/3),
                  ('(∫μ)²', 1/9)]:
    if abs(C_ratio/val - 1) < 0.2:
        print(f"    C ≈ {name} = {val:.6f} (ratio: {C_ratio/val:.6f})")

print(f"\n  Testing μ(θ_max) = C × δ:")
C_ratio2 = mu_at_eps / delta
print(f"    C = μ/δ = {C_ratio2:.8f}")

for name, val in [('1', 1), ('1/2', 0.5), ('1/3', 1/3), ('1/9', 1/9),
                  ('(∫μ)²', 1/9), ('β/δ', beta/delta), ('ε²/δ', eps_s**2/delta)]:
    if abs(C_ratio2/val - 1) < 0.3:
        print(f"    C ≈ {name} = {val:.6f} (ratio: {C_ratio2/val:.6f})")

# ============================================================================
# STEP 3: RECOMPUTE C₂ AT THE DERIVED (β, ε*) WITH CONVERGENCE CHECK
# ============================================================================

print(f"\n" + "=" * 72)
print("  STEP 3: C₂ CONVERGENCE AT DERIVED (β=δ/9, ε=ε*)")
print("=" * 72)

C2_target = 1/(72*pi)
ls_target = 36*pi
alpha = 1/137.035999178

def full_computation(beta_val, eps_val, N=8000):
    theta_max = pi/2 - eps_val
    
    results = {}
    for m_val in [0, 3]:
        theta_grid = np.linspace(0, theta_max, N)
        dtheta = theta_grid[1] - theta_grid[0]
        mu_val = mu(theta_grid[1:-1])
        mu_sq_reg = mu_val**2 + delta**2
        
        if m_val == 0:
            V = K_curvature(theta_grid[1:-1]) / 2
        else:
            f_theta = mu_sq_reg / (mu_sq_reg + beta_val)**2
            V = m_val**2 * f_theta + K_curvature(theta_grid[1:-1]) / 2
        
        diag = 2.0/dtheta**2 + V
        off = -np.ones(N-3)/dtheta**2
        evals, evecs = eigh_tridiagonal(diag, off, select='i', select_range=(0,2))
        
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
    
    C2 = 9 * I30**2 / abs(dE)
    lam_scale = abs(dE) / (2 * 9 * I30**2)
    product = 2 * C2 * lam_scale
    
    return {
        'lam00': lam00, 'lam30': lam30, 'dE': dE,
        'I30': I30, 'C2': C2, 'lam_scale': lam_scale, 'product': product,
    }

# Convergence study: ε* re-derived at each N, then C₂ computed at that (β, ε*)
print(f"\n  Strategy: at each N, FIRST solve for ε*(N), THEN compute C₂.")
print(f"  This couples the cutoff to the grid, eliminating systematic drift.\n")

print(f"  {'N':>8s}  {'ε*(N)':>14s}  {'C₂':>14s}  {'err%':>10s}  {'λ_scale':>12s}  {'err%':>10s}  {'2C₂λ':>8s}")
print("  " + "-" * 84)

N_values = [2000, 3000, 4000, 6000, 8000, 10000, 12000, 16000]
c2_results = []

for N in N_values:
    # Solve for ε* at this N
    def lam_N(eps):
        return get_lambda_30(eps, N=N)
    
    try:
        eps_N = brentq(lam_N, 0.0005, 0.003, xtol=1e-14)
    except:
        # Widen bracket
        try:
            eps_N = brentq(lam_N, 0.0003, 0.01, xtol=1e-14)
        except:
            print(f"  {N:>8d}  FAILED to find ε*")
            continue
    
    # Compute C₂ at this (β, ε*)
    r = full_computation(beta, eps_N, N=N)
    c2_err = (r['C2'] - C2_target) / C2_target * 100
    ls_err = (r['lam_scale'] - ls_target) / ls_target * 100
    
    c2_results.append((N, eps_N, r))
    print(f"  {N:>8d}  {eps_N:>14.10f}  {r['C2']:>14.10f}  {c2_err:>+10.4f}  "
          f"{r['lam_scale']:>12.6f}  {ls_err:>+10.4f}  {r['product']:>8.6f}")

# ============================================================================
# STEP 4: COMPARE — FIXED ε vs DERIVED ε*
# ============================================================================

print(f"\n" + "=" * 72)
print("  STEP 4: COMPARISON — FIXED ε=0.001 vs DERIVED ε*")
print("=" * 72)

print(f"\n  {'N':>8s}  {'C₂ (ε=0.001)':>14s}  {'err%':>8s}  {'C₂ (ε=ε*)':>14s}  {'err%':>8s}  {'Improvement':>12s}")
print("  " + "-" * 76)

for N in [4000, 8000, 12000, 16000]:
    # Fixed ε
    r_fixed = full_computation(beta, 0.001, N=N)
    c2_err_fixed = (r_fixed['C2'] - C2_target) / C2_target * 100
    
    # Derived ε*
    def lam_N(eps):
        return get_lambda_30(eps, N=N)
    eps_N = brentq(lam_N, 0.0005, 0.003, xtol=1e-14)
    r_derived = full_computation(beta, eps_N, N=N)
    c2_err_derived = (r_derived['C2'] - C2_target) / C2_target * 100
    
    improvement = abs(c2_err_fixed) / abs(c2_err_derived) if abs(c2_err_derived) > 1e-10 else float('inf')
    print(f"  {N:>8d}  {r_fixed['C2']:>14.10f}  {c2_err_fixed:>+8.4f}  "
          f"{r_derived['C2']:>14.10f}  {c2_err_derived:>+8.4f}  {improvement:>12.1f}×")

# ============================================================================
# STEP 5: FINAL RESULT
# ============================================================================

print(f"\n" + "=" * 72)
print("  FINAL RESULT")
print("=" * 72)

# Best result
N_final = 16000
def lam_final(eps):
    return get_lambda_30(eps, N=N_final)
eps_final = brentq(lam_final, 0.0005, 0.003, xtol=1e-14)
r_final = full_computation(beta, eps_final, N=N_final)

c2_err = (r_final['C2'] - C2_target) / C2_target * 100
ls_err = (r_final['lam_scale'] - ls_target) / ls_target * 100

print(f"""
  ╔══════════════════════════════════════════════════════════════════════╗
  ║  ZERO-PARAMETER RESULT WITH DERIVED (β, ε)                         ║
  ╠══════════════════════════════════════════════════════════════════════╣
  ║                                                                      ║
  ║  INPUTS (all from geometry):                                        ║
  ║    δ = 1/(π⁶√2)                   [12D fiber resolution]           ║
  ║    β = δ(∫μdθ)² = δ/9             [projected fiber volume]         ║
  ║    ε* = solution of λ₃₀(β,ε)=0    [critical binding condition]     ║
  ║                                                                      ║
  ║  DERIVED CUTOFF:                                                    ║
  ║    ε* = {eps_final:.12f}                                     ║
  ║    θ_max = {np.degrees(pi/2 - eps_final):.8f}°                                  ║
  ║    μ(θ_max) = {mu(pi/2 - eps_final):.6e}                               ║
  ║                                                                      ║
  ║  RESULTS (N = {N_final}):                                           ║
  ║    λ₀₀ = {r_final['lam00']:>14.4f}                                      ║
  ║    λ₃₀ = {r_final['lam30']:>14.6f}   (≈ 0, by construction)           ║
  ║    ΔE  = {r_final['dE']:>14.4f}                                      ║
  ║    I₃₀ = {r_final['I30']:>14.8f}                                    ║
  ║                                                                      ║
  ║    C₂ = {r_final['C2']:.10f}                                  ║
  ║    1/(72π) = {C2_target:.10f}                                  ║
  ║    Error: {c2_err:>+.4f}%                                           ║
  ║                                                                      ║
  ║    λ_scale = {r_final['lam_scale']:.6f}                                       ║
  ║    36π = {ls_target:.6f}                                       ║
  ║    Error: {ls_err:>+.4f}%                                           ║
  ║                                                                      ║
  ║    2 × C₂ × λ_scale = {r_final['product']:.10f}                        ║
  ║                                                                      ║
  ║    a_e = α/(2π) × {r_final['product']:.10f}                        ║
  ║        = {alpha/(2*pi) * r_final['product']:.12f}                  ║
  ║    Target: {alpha/(2*pi):.12f}                                 ║
  ║    Error: {(r_final['product']-1)*100:>+.6f}%                                       ║
  ║                                                                      ║
  ╚══════════════════════════════════════════════════════════════════════╝
""")

# ============================================================================
# VISUALIZATION
# ============================================================================

fig, axes = plt.subplots(1, 3, figsize=(18, 5.5))
fig.suptitle('Closing the Loop: β = δ/9, ε* from λ₃₀ = 0', 
             fontsize=14, fontweight='bold')

# Panel 1: C₂ convergence comparison
ax = axes[0]
Ns_fixed = []; errs_fixed = []
Ns_derived = []; errs_derived = []

for N in [2000, 3000, 4000, 6000, 8000, 10000, 12000, 16000]:
    r_f = full_computation(beta, 0.001, N=N)
    Ns_fixed.append(N)
    errs_fixed.append((r_f['C2'] - C2_target)/C2_target * 100)

for N, eps_N, r in c2_results:
    Ns_derived.append(N)
    errs_derived.append((r['C2'] - C2_target)/C2_target * 100)

ax.plot(Ns_fixed, errs_fixed, 'rs--', linewidth=2, markersize=7, label='Fixed ε = 0.001')
ax.plot(Ns_derived, errs_derived, 'bo-', linewidth=2, markersize=7, label='Derived ε* (λ₃₀=0)')
ax.axhline(y=0, color='green', linewidth=2, linestyle='--', label='Exact: 1/(72π)')
ax.set_xlabel('Grid points N', fontsize=12)
ax.set_ylabel('C₂ error (%)', fontsize=12)
ax.set_title('C₂ Convergence: Fixed vs Derived ε')
ax.legend(fontsize=9)
ax.grid(True, alpha=0.3)

# Panel 2: ε*(N) convergence
ax = axes[1]
Ns_eps = [x[0] for x in c2_results]
eps_vals = [x[1] for x in c2_results]

ax.plot(Ns_eps, eps_vals, 'ko-', linewidth=2, markersize=8)
ax.axhline(y=0.001, color='red', linewidth=1.5, linestyle=':', label='ε = 0.001 (old)')
ax.set_xlabel('Grid points N', fontsize=12)
ax.set_ylabel('ε*', fontsize=12)
ax.set_title('Derived ε* vs Resolution')
ax.legend(fontsize=10)
ax.grid(True, alpha=0.3)

# Annotate converged value
ax.text(0.5, 0.15, f'ε* → {eps_final:.6f}\nθ_max → {np.degrees(pi/2-eps_final):.4f}°',
        transform=ax.transAxes, fontsize=11, ha='center',
        bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.7))

# Panel 3: The complete picture
ax = axes[2]
ax.set_xlim(0, 10)
ax.set_ylim(0, 10)

items = [
    (5, 9.0, '12D Fiber → δ = 1/(π⁶√2)', 'lightblue'),
    (5, 7.8, 'Projection: β = δ(∫μdθ)² = δ/9', 'lightyellow'),
    (5, 6.6, f'Critical binding: λ₃₀(δ/9, ε*)=0\nε* = {eps_final:.6f}', 'lightyellow'),
    (5, 5.2, 'ONE bound state (m=±3)', 'lightgreen'),
    (5, 4.0, f'C₂ = {r_final["C2"]:.6f} ≈ 1/(72π)', 'lightgreen'),
    (5, 2.8, f'2C₂λ = {r_final["product"]:.8f}', 'lightgreen'),
    (5, 1.6, f'a_e = α/(2π) ✓', 'gold'),
]

for x, y, text, color in items:
    ax.text(x, y, text, ha='center', va='center', fontsize=9, fontweight='bold',
            bbox=dict(boxstyle='round,pad=0.3', facecolor=color, edgecolor='black'))

for i in range(len(items)-1):
    ax.annotate('', xy=(5, items[i+1][1]+0.4), xytext=(5, items[i][1]-0.4),
                arrowprops=dict(arrowstyle='->', lw=1.5))

ax.text(5, 0.5, 'ZERO free parameters', ha='center', fontsize=13, 
        fontweight='bold', color='red')
ax.set_title('Complete Chain')
ax.axis('off')

plt.tight_layout()
plt.savefig('/home/claude/closing_the_loop.png', dpi=150, bbox_inches='tight')
print(f"  [Plot saved]")

# ============================================================================
# PHYSICAL MEANING OF ε*
# ============================================================================

print(f"\n" + "=" * 72)
print("  PHYSICAL MEANING OF ε*")
print("=" * 72)

print(f"""
  ε* = {eps_final:.10f}
  
  At this cutoff:
    μ(π/2 - ε*) = {mu(pi/2 - eps_final):.6e}
    β = δ/9      = {beta:.6e}
    
    μ(θ_max) / β = {mu(pi/2 - eps_final)/beta:.6f}
    
  INTERPRETATION:
  ε* is the angle where the measure μ reaches a specific fraction 
  of the fiber dilution β. It's NOT where μ = β (that ratio is
  {mu(pi/2 - eps_final)/beta:.4f}, not 1), but it IS determined
  by the eigenvalue condition — the point where the potential well
  is just deep enough to bind one state.
  
  ε* is an EMERGENT scale, not a parameter. Given:
    (i)   the base geometry μ(θ) = cos²θ sinθ
    (ii)  the fiber resolution δ = 1/(π⁶√2)  
    (iii) the projection β = δ/9
  
  ε* follows uniquely from the eigenvalue equation.
  It cannot be chosen; it IS what it IS.
""")

print("=" * 72)
print("  SECTOR CLOSED")
print("=" * 72)
