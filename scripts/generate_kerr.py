# scripts/generate_kerr.py
import sympy as sp
from sympy.printing.cxx import cxxcode

print("\n//KerrChristoffel.inl")
# coordinates and constants
t, x, y, z = sp.symbols('t x y z')
coords = [t, x, y, z]
M, a = sp.symbols('M a')

R2 = x**2 + y**2 + z**2
a2 = a**2

half_term = (R2 - a2) * 0.5
r2 = half_term + sp.sqrt(half_term**2 + a2 * z**2)
r = sp.sqrt(r2)

f = (2.0 * M * r2 * r) / (r2**2 + a2 * z**2)

# Covariant null vector l_mu (Index DOWN)
l_down = sp.Matrix([
    [1.0],
    [(r * x + a * y) / (r2 + a2)],
    [(r * y - a * x) / (r2 + a2)],
    [z / r]
])

# flat space (eta)
eta = sp.diag(-1.0, 1.0, 1.0, 1.0)
eta_inv = sp.diag(-1.0, 1.0, 1.0, 1.0) # Inverse of eta is itself

g = eta + f * (l_down * l_down.T) # g_mu_nu
l_up = eta_inv * l_down

# FAST EXACT INVERSE: g^{mu nu} = eta^{mu nu} - f * l^mu * l^nu
g_inv = eta_inv - f * (l_up * l_up.T)

print("Calculating derivatives...")
gamma_list = []
gamma_symbols = []

for mu in range(4):
    for alpha in range(4):
        for beta in range(4):
            sum_term = 0
            for lam in range(4):
                term = sp.diff(g[beta, lam], coords[alpha]) + \
                       sp.diff(g[alpha, lam], coords[beta]) - \
                       sp.diff(g[alpha, beta], coords[lam])
                sum_term += g_inv[mu, lam] * term
            
            expr = 0.5 * sum_term
            gamma_list.append(expr)
            gamma_symbols.append(f"gamma[{mu}]({alpha}, {beta})")

print("Optimising into C++...")
replacements, reduced_exprs = sp.cse(gamma_list)

print("\n// --- \n")

# temporary variables (x0 = ..., x1 = ...)
for var, expr in replacements:
    print(f"double {var} = {cxxcode(expr)};")

print("\n// Christoffel Assignments:")
# FINAL assignments
for i, expr in enumerate(reduced_exprs):
    if expr != 0:
        print(f"{gamma_symbols[i]} = {cxxcode(expr)};")

print("\n// --- ")
print("\n//KerrMetric.inl")

g_list = []
g_symbols = []

for mu in range(4):
    for nu in range(4):
        g_list.append(g[mu, nu])
        g_symbols.append(f"g({mu}, {nu})")

print("Optimising Metric into C++...")
replacements_g, reduced_exprs_g = sp.cse(g_list)

# Print temporary variables for the metric
for var, expr in replacements_g:
    print(f"double {var} = {cxxcode(expr)};")

print("\n// Metric Assignments:")
for i, expr in enumerate(reduced_exprs_g):
    if expr != 0:
        print(f"{g_symbols[i]} = {cxxcode(expr)};")

print("\n//KerrShadow.inl")
print("// Analytical limits for the photon region and Bardeen impact parameters")

# The constants of motion xi (Lz/E) and eta (Q/E^2) parameterised by photon orbit radius r_p
r_p = sp.Symbol('r_p')
xi = (r_p**2 * (r_p - 3*M) + a**2 * (r_p + M)) / (a * (r_p - M))
eta = (r_p**3 * (4*M*a**2 - r_p * (r_p - 3*M)**2)) / (a**2 * (r_p - M)**2)

# Generate C++ code for xi and eta functions
print("inline double compute_xi(double r_p, double M, double a) {")
print(f"    return {cxxcode(xi)};")
print("}")

print("inline double compute_eta(double r_p, double M, double a) {")
print(f"    return {cxxcode(eta)};")
print("}")

