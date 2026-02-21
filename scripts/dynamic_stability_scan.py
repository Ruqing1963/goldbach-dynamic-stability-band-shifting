"""
dynamic_stability_scan.py — Full scan [1024, 2048] for Paper #10
Tracks: ρ evolution, band shifting law, non-emptiness, rad_odd(M) correlation
"""
import math
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import csv

def sieve_and_radicals(limit):
    is_prime = bytearray([1]) * (limit + 1)
    is_prime[0] = is_prime[1] = 0
    rads = np.ones(limit + 1, dtype=np.int64)
    for i in range(2, limit + 1):
        if is_prime[i]:
            for j in range(i, limit + 1, i):
                if j != i: is_prime[j] = 0
                rads[j] *= i
    return is_prime, rads

def odd_rad(x, rads):
    r = int(rads[x])
    while r % 2 == 0: r //= 2
    return r if r > 0 else 1

LIMIT = 20000
is_prime, rads = sieve_and_radicals(LIMIT)

# ═══════════════════════════════════════════════════════════════════════════════
# FULL SCAN: every even N in [1024, 2048]
# ═══════════════════════════════════════════════════════════════════════════════
results = []
for N in range(1024, 2050, 2):
    M = N // 2
    rM = odd_rad(M, rads)
    log_rM = math.log(rM) if rM > 1 else 0.0

    gb_rhos, comp_rhos, mix_rhos = [], [], []
    gb_ground_p = None

    for p in range(3, N // 2 + 1):
        q = N - p
        if q <= 1: continue
        rp = odd_rad(p, rads)
        rq = odd_rad(q, rads)
        base = rp * rq * rM * rM
        cond = base * base
        rho = math.log(cond) / math.log(N) if cond > 1 else 0.0

        if is_prime[p] and is_prime[q]:
            gb_rhos.append(rho)
            if gb_ground_p is None or rho < min(gb_rhos[:-1], default=999):
                gb_ground_p = p
        elif not is_prime[p] and not is_prime[q]:
            comp_rhos.append(rho)
        else:
            mix_rhos.append(rho)

    results.append({
        'N': N, 'M': M, 'rad_M': rM, 'log_rad_M': log_rM,
        'num_gb': len(gb_rhos),
        'rho_min': min(gb_rhos) if gb_rhos else None,
        'rho_mean': np.mean(gb_rhos) if gb_rhos else None,
        'rho_max': max(gb_rhos) if gb_rhos else None,
        'bw_gb': (max(gb_rhos) - min(gb_rhos)) if len(gb_rhos) >= 2 else 0,
        'comp_rhos': comp_rhos,
        'ground_p': gb_ground_p,
    })

# ═══════════════════════════════════════════════════════════════════════════════
# ANALYSIS 1: Non-emptiness
# ═══════════════════════════════════════════════════════════════════════════════
total = len(results)
nonempty = sum(1 for r in results if r['num_gb'] > 0)
print("=" * 80)
print("ANALYSIS 1: Non-emptiness in [1024, 2048]")
print("=" * 80)
print(f"  Total even N scanned: {total}")
print(f"  N with ≥ 1 Goldbach pair: {nonempty}")
print(f"  N with 0 Goldbach pairs: {total - nonempty}")
print(f"  Non-emptiness rate: {nonempty/total*100:.1f}%")

# ═══════════════════════════════════════════════════════════════════════════════
# ANALYSIS 2: Band Shifting Law — ρ_mean vs log(rad_odd(M))
# ═══════════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 80)
print("ANALYSIS 2: Band Shifting Law")
print("=" * 80)

valid = [r for r in results if r['rho_mean'] is not None and r['log_rad_M'] > 0]
log_rMs = np.array([r['log_rad_M'] for r in valid])
rho_means = np.array([r['rho_mean'] for r in valid])
rho_mins = np.array([r['rho_min'] for r in valid])

# Linear regression: rho_mean = a * log(rad_M) / log(N) + b
# Normalise: define the "conduit contribution" as 2*log(rad_M)/log(N)
conduit_contrib = np.array([2 * r['log_rad_M'] / math.log(r['N']) for r in valid])
# Predict: rho_mean ≈ conduit_contrib + boundary_contribution
boundary = rho_means - conduit_contrib

from numpy.polynomial import polynomial as P
# Fit rho_mean vs conduit_contrib
coeffs = np.polyfit(conduit_contrib, rho_means, 1)
slope, intercept = coeffs
residuals = rho_means - (slope * conduit_contrib + intercept)
r_squared = 1 - np.var(residuals) / np.var(rho_means)

print(f"  Linear fit: ⟨ρ⟩_GB = {slope:.4f} × [2log(rad_M)/log(N)] + {intercept:.4f}")
print(f"  R² = {r_squared:.6f}")
print(f"  Residual std = {np.std(residuals):.4f}")

# Fit rho_min vs conduit_contrib
coeffs_min = np.polyfit(conduit_contrib, rho_mins, 1)
slope_min, intercept_min = coeffs_min
residuals_min = rho_mins - (slope_min * conduit_contrib + intercept_min)
r2_min = 1 - np.var(residuals_min) / np.var(rho_mins)
print(f"\n  Linear fit: ρ_min = {slope_min:.4f} × [2log(rad_M)/log(N)] + {intercept_min:.4f}")
print(f"  R² = {r2_min:.6f}")

# ═══════════════════════════════════════════════════════════════════════════════
# ANALYSIS 3: Jump statistics
# ═══════════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 80)
print("ANALYSIS 3: Step-to-step jump statistics")
print("=" * 80)

jumps_rho_min = []
jumps_conduit = []
for i in range(1, len(results)):
    if results[i]['rho_min'] is not None and results[i-1]['rho_min'] is not None:
        jumps_rho_min.append(results[i]['rho_min'] - results[i-1]['rho_min'])
        c_now = 2 * results[i]['log_rad_M'] / math.log(results[i]['N'])
        c_prev = 2 * results[i-1]['log_rad_M'] / math.log(results[i-1]['N'])
        jumps_conduit.append(c_now - c_prev)

jumps_rho_min = np.array(jumps_rho_min)
jumps_conduit = np.array(jumps_conduit)

print(f"  Δρ_min per step: mean={np.mean(jumps_rho_min):.4f}, std={np.std(jumps_rho_min):.4f}")
print(f"  Δρ_min range: [{min(jumps_rho_min):.4f}, {max(jumps_rho_min):.4f}]")
print(f"  |Δρ_min| > 1: {np.sum(np.abs(jumps_rho_min) > 1)} out of {len(jumps_rho_min)} steps")
print(f"  |Δρ_min| > 2: {np.sum(np.abs(jumps_rho_min) > 2)} out of {len(jumps_rho_min)} steps")

# Correlation between jump in ρ and jump in conduit
corr = np.corrcoef(jumps_conduit, jumps_rho_min)[0, 1]
print(f"\n  Correlation(Δconduit, Δρ_min) = {corr:.4f}")

# ═══════════════════════════════════════════════════════════════════════════════
# ANALYSIS 4: Bandwidth stability
# ═══════════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 80)
print("ANALYSIS 4: Goldbach bandwidth statistics")
print("=" * 80)
bws = [r['bw_gb'] for r in results if r['num_gb'] >= 2]
print(f"  BW_GB: mean={np.mean(bws):.4f}, std={np.std(bws):.4f}")
print(f"  BW_GB range: [{min(bws):.4f}, {max(bws):.4f}]")

# ═══════════════════════════════════════════════════════════════════════════════
# FIGURE 1: Propagation overview [1024, 2048]
# ═══════════════════════════════════════════════════════════════════════════════
print("\n[*] Generating figures...")

fig = plt.figure(figsize=(14, 10))
gs = GridSpec(3, 1, height_ratios=[3, 1, 1.2], hspace=0.08)

ax1 = fig.add_subplot(gs[0])
ax2 = fig.add_subplot(gs[1], sharex=ax1)
ax3 = fig.add_subplot(gs[2], sharex=ax1)

Ns = [r['N'] for r in results if r['rho_mean'] is not None]
rmins = [r['rho_min'] for r in results if r['rho_mean'] is not None]
rmeans_p = [r['rho_mean'] for r in results if r['rho_mean'] is not None]
rmaxs = [r['rho_max'] for r in results if r['rho_mean'] is not None]
radMs = [r['rad_M'] for r in results if r['rho_mean'] is not None]
nGBs = [r['num_gb'] for r in results if r['rho_mean'] is not None]

# Top: ρ band
ax1.fill_between(Ns, rmins, rmaxs, alpha=0.12, color='#2255BB')
ax1.plot(Ns, rmeans_p, '-', color='#2255BB', linewidth=0.6, alpha=0.7,
         label='$\\langle\\rho\\rangle_{\\mathrm{GB}}$')
ax1.plot(Ns, rmins, '-', color='#CC3333', linewidth=0.6, alpha=0.7,
         label='$\\rho_{\\min}$')

# Mark 2^k points
for k in [10, 11]:
    Nk = 2**k
    idx = next(i for i, r in enumerate(results) if r['N'] == Nk)
    if results[idx]['rho_min'] is not None:
        ax1.plot(Nk, results[idx]['rho_min'], 'v', color='red', markersize=10, zorder=5)
        ax1.annotate(f'$2^{{{k}}}$', xy=(Nk, results[idx]['rho_min']),
                     xytext=(Nk + 15, results[idx]['rho_min'] - 0.3),
                     fontsize=10, color='red', fontweight='bold')

ax1.set_ylabel("Chen's Ratio $\\rho$", fontsize=12)
ax1.set_title('Conductor Orbit Propagation: $N \\in [2^{10},\\, 2^{11}]$', fontsize=14)
ax1.legend(loc='lower right', fontsize=9)
ax1.grid(True, alpha=0.12)
plt.setp(ax1.get_xticklabels(), visible=False)

# Middle: rad_odd(M)
ax2.bar(Ns, radMs, width=1.5, color='#CC6644', alpha=0.5)
ax2.set_ylabel('$\\mathrm{rad}_{\\mathrm{odd}}(M)$', fontsize=10)
ax2.set_yscale('log')
ax2.grid(True, alpha=0.12)
plt.setp(ax2.get_xticklabels(), visible=False)

# Bottom: #GB pairs
ax3.bar(Ns, nGBs, width=1.5, color='#228833', alpha=0.5)
ax3.set_ylabel('\\# Goldbach pairs', fontsize=10)
ax3.set_xlabel('$N$', fontsize=12)
ax3.grid(True, alpha=0.12)

plt.savefig('/home/claude/fig1_propagation_full.pdf', dpi=300, bbox_inches='tight')
plt.savefig('/home/claude/fig1_propagation_full.png', dpi=200, bbox_inches='tight')
plt.close()
print("  Saved fig1_propagation_full.pdf/.png")

# ═══════════════════════════════════════════════════════════════════════════════
# FIGURE 2: Band Shifting Law — ρ_mean vs conduit contribution
# ═══════════════════════════════════════════════════════════════════════════════
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5.5))

# Left: scatter ρ_mean vs 2log(rad_M)/log(N)
cc = np.array([2 * r['log_rad_M'] / math.log(r['N']) for r in valid])
ax1.scatter(cc, rho_means, s=6, alpha=0.4, color='#2255BB', label='$\\langle\\rho\\rangle_{\\mathrm{GB}}$')
ax1.scatter(cc, rho_mins, s=4, alpha=0.3, color='#CC3333', label='$\\rho_{\\min}$')

# Regression line
xx = np.linspace(0, max(cc) * 1.05, 100)
ax1.plot(xx, slope * xx + intercept, 'k--', linewidth=1.2,
         label=f'Fit: $\\langle\\rho\\rangle = {slope:.2f}\\,\\xi + {intercept:.2f}$\n$R^2 = {r_squared:.4f}$')

ax1.set_xlabel('Static conduit contribution $\\xi = 2\\log(\\mathrm{rad}_{\\mathrm{odd}}(M))/\\log N$', fontsize=10)
ax1.set_ylabel("Chen's Ratio $\\rho$", fontsize=11)
ax1.set_title('Band Shifting Law', fontsize=13)
ax1.legend(fontsize=9)
ax1.grid(True, alpha=0.15)

# Right: residual histogram
ax2.hist(residuals, bins=40, color='#2255BB', alpha=0.6, edgecolor='white')
ax2.axvline(x=0, color='black', linestyle='--', linewidth=0.8)
ax2.set_xlabel('Residual $\\langle\\rho\\rangle - \\hat{\\rho}$', fontsize=11)
ax2.set_ylabel('Count', fontsize=11)
ax2.set_title(f'Residuals (std = {np.std(residuals):.3f})', fontsize=12)
ax2.grid(True, alpha=0.15)

plt.suptitle('Figure 2 — Band Shifting Law: $\\langle\\rho\\rangle_{\\mathrm{GB}}$ vs. Static Conduit Contribution',
             fontsize=13, y=1.02)
plt.tight_layout()
plt.savefig('/home/claude/fig2_band_shifting.pdf', dpi=300, bbox_inches='tight')
plt.savefig('/home/claude/fig2_band_shifting.png', dpi=200, bbox_inches='tight')
plt.close()
print("  Saved fig2_band_shifting.pdf/.png")

# ═══════════════════════════════════════════════════════════════════════════════
# Export CSV
# ═══════════════════════════════════════════════════════════════════════════════
with open('/home/claude/propagation_data.csv', 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['N', 'rad_odd_M', 'num_goldbach', 'rho_min', 'rho_mean', 'rho_max', 'BW_goldbach'])
    for r in results:
        if r['rho_mean'] is not None:
            w.writerow([r['N'], r['rad_M'], r['num_gb'],
                        f"{r['rho_min']:.6f}", f"{r['rho_mean']:.6f}",
                        f"{r['rho_max']:.6f}", f"{r['bw_gb']:.6f}"])
print("  Saved propagation_data.csv")
print("\n[*] Done.")
