# ------------------------------------------------------------------------------
# ---- 1. INTRODUCTION AND CONTENTS --------------------------------------------

# The last time I rewrite the code. Promise. 

#= Use this line in the Terminal to start Julia. 
/Applications/Julia-1.8.app/Contents/Resources/julia/bin/julia
=#

# Packages used in this code.
using Dates                         
using LinearAlgebra
using SpecialFunctions              # 2.4.0 
using Random
using IntervalArithmetic            # 0.22.14
using QuadGK                        # 2.9.4
using Readables                     # 0.3.3
using IntervalArithmetic.Symbols

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------



# ------------------------------------------------------------------------------
# ---- 2. BASIC FORMULAS -------------------------------------------------------

# BASIC FORMULAS

I(x, y) = interval(@interval(x), @interval(y))

# "@interval" is back and works as it should: 
# 
# julia> BigFloat.(bounds(I(0.1)))
# (0.09999999999999999167332731531132594682276248931884765625, 0.1000000000000000055511151231257827021181583404541015625)

# The "squaring" function. 

# IntervalArithmetic uses BigFloat numbers with powers of intervals, so computing "X^2" is much slower than "X * X". 
# "Note, in particular, that in order to obtain correct rounding for the power function (^), intervals are converted to and from BigFloat; this implies a significant slow-down in this case."
# https://juliaintervals.github.io/IntervalArithmetic.jl/latest/usage/
# Since we square many quantities in our functions below, we define a function "SQ" which is defined by SQ(x) = x * x. 

SQ(X) = @interval X * X

function SQI(X) 
    l = inf(X)
    r = sup(X)
    if r < 0
        return I(r^2, l^2)
    elseif l > 0
        return I(l^2, r^2)
    else 
        return I(0, max(l, r)^2)
    end 
end

# The volume of the d-dimensional unit ball.
    K(d) = @interval(π^(d / 2)) / (gamma(d / 2 + 1) ± eps())
κ_big(d) = big(π)^(d / 2) / gamma(big(d / 2 + 1))

# The function B (formerly Υ).
B(A) = @interval acos((1 / 4) * SQ(cos(A)) + cos(A) * √((1 / 16) * SQ(cos(A)) + 1 / 2))

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------



# ------------------------------------------------------------------------------
# ---- 3. THRESHOLD ANGLES; A and A★ -------------------------------------------

A_test = [[ I(0, 0)     I(0, 0)     ]
          [ I(1.3, 1.4) I(0, 0)     ]
          [ I(1.2, 1.3) I(3.3, 3.4) ]]

# Converts a vector A to a triangular matrix. 
function ◺(A)
    if length(A) == 6 # If the configuration points are in ℝ^4, then A has six entries.
        return [[0    0    0   ]
                [A[1] 0    0   ]
                [A[2] A[3] 0   ]
                [A[4] A[5] A[6]]]
    else # If the configuration points are in ℝ^3, then A has three entries.
        return [[0    0   ]
                [A[1] 0   ]
                [A[2] A[3]]]
    end
end

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------



# ------------------------------------------------------------------------------
# ---- 4. THE DISTANCE FROM THE ORIGIN TO A FACE; H AND H̃ ---------------------

# The lower bounds for H.

H_min_2(A★) = min(IntervalArithmetic.interval(2 * sin(A★[2] / 2)), 
                   IntervalArithmetic.interval(2 * sin(A★[3] / 2)), IntervalArithmetic.interval(2 * sin(A★[4] / 2)))

# Lower bound for the distance from the origin to the plane aff{2y^i, 2y^j, 2y^k}. 

function H123(A)
    C1 = @interval 1
    C2 = @interval (1 - cos(A[2,1])) / sin(A[2,1])
    C3 = @interval 1 / sin(A[3,1]) + 1 / tan(A[3,1]) * (sin(A[3,2] - A[2,1]) - sin(A[3,2])) / sin(A[2,1])
    return @interval 2 / √(SQ(C1) + SQ(C2) + SQ(C3))
end

# The lower bounds for H.

function H123_min(A★)
    C1 = @interval 1
    C2 = @interval (1 + cos(B(A★[3]))) / sin(B(A★[3]))
    C3 = @interval 1 / sin(A★[3]) + 1 / (sin(B(A★[3]) / 2) * tan(A★[3]))
    return @interval 2 / √(SQ(C1) + SQ(C2) + SQ(C3))
end

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------



# ------------------------------------------------------------------------------
# ---- 5. THE CONFIGURATION POINTS Y^i IN TERMS OF THE ANGLES A ----------------

# Notation:
#   "Yij" = Y^j in ℝ^i. Note that Yij is an interval. 

# The points Y^1, Y^2, and Y^3 in three dimensions.
Y31(A) = [ @interval(1), 
           @interval(0), 
           @interval(0) ]
Y32(A) = [ @interval(cos(A[2,1])),
           @interval(sin(A[2,1])),
           @interval(0) ]
Y33(A) = [ @interval(cos(A[3,1]) * cos(A[3,2])), 
           @interval(cos(A[3,1]) * sin(A[3,2])), 
           @interval(sin(A[3,1])) ] 

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------



# ------------------------------------------------------------------------------
# ---- 6. THE OUTWARD NORMAL VECTORS; Uijk -------------------------------------

# The outward normal vectors to the faces F_ijk, calculated from the null spaces of {Y^i, Y^j, Y^k}.

# Notation:
#   "Nijk" = the outward normal to the face F_ijk of P^4.
#   "Uijk" = the normalized outward normal to the face F_ijk of P^4.
#   "_n" = normalizing nijk using the norm function.
#   "_t" = normalizing nijk using trigonometric functions.

# Note that the _t expressions result in narrower intervals than the _n expressions.

N123(A) = [ IntervalArithmetic.interval(  0 ),  
            IntervalArithmetic.interval(  0 ),  
            IntervalArithmetic.interval(  0 ),            
            IntervalArithmetic.interval( -1 ) ] 
            # FFS I had "1" instead of "-1" and spent a day debugging this sh*t 

# The normalized outward normal vectors using the "norm" function.
U123_n(A) = N123(A)#/ norm(N123(A)) 

# The denominator of Uijk — excluding the square root — expressed in trigonometric functions.
D123(A) = IntervalArithmetic.interval(1)

# The normalized outward normal vectors using trigonometric functions.
U123_t(A) = N123(A)#/ √(D123(A))

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------



# ------------------------------------------------------------------------------
# ---- 7. THE ARC LENGTH; fij AND gij ------------------------------------------

# Notation:
#     "Fij_m" = The angle F_ij between the vertices Y^i and Y^j 
#               in m dimensions.
#     "Gkl_m" = The angle G_kl between the points U_ijk and U_ijl 
#               in m dimensions.
#     "_nd"   = "norm and dot," these expressions are written 
#               as an inner product of the terms uijk_n and uijl_n.
#     "_t"    = "trigonometric," these expressions are written 
#               in terms of trigonometric functions only.

# The _t expressions result in narrower intervals than the _nd expressions. 

# The angles F_ij between the vertices Y^i and Y^j in three dimensions.
F12_3_d(A) = @interval(A[2,1])
F13_3_d(A) = @interval(acos(Y31(A) ⋅ Y33(A)))
F23_3_d(A) = @interval(acos(Y32(A) ⋅ Y33(A)))
F12_3_t(A) = @interval(A[2,1])
F13_3_t(A) = @interval(acos(cos(A[3,1]) * cos(A[3,2])))
F23_3_t(A) = @interval(acos(cos(A[3,1]) * cos(A[3,2] - A[2,1])))

# The angles G_jk between the normal vectors U_ijk and U_ijl in three dimensions.
G12_3_t(A) = @interval(acos((-SQ(sin(A[3,1])) * cos(A[2,1]) - SQ(cos(A[3,1])) * sin(A[3,2]) * sin(A[3,2] - A[2,1])) / (√(1 - SQ(cos(A[3,1])) * SQ(cos(A[3,2]))) * √(1 - SQ(cos(A[3,1])) * SQ(cos(A[3,2] - A[2,1]))))))
G13_3_t(A) = @interval(acos((cos(A[3,1]) * sin(A[3,2] - A[2,1])) / √(1 - SQ(cos(A[3,1])) * SQ(cos(A[3,2] - A[2,1])))))
G23_3_t(A) = @interval(acos((-cos(A[3,1]) * sin(A[3,2])) / √(1 - SQ(cos(A[3,1])) * SQ(cos(A[3,2])))))

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------



# ------------------------------------------------------------------------------
# ---- 8. THE AREA OF A SPHERICAL TRIANGLE; fijk AND gjkl ----------------------

# The solid angle f_123 formed by the vertices y^1, y^2, and y^3 in three dimensions. This formula does not use L'Huilier's theorem. 
F123_3(A) = @interval(2 * π) - (G23_3_t(A) + G13_3_t(A) + G12_3_t(A))
G123_3(A) = @interval(2 * π) - (F12_3_t(A) + F13_3_t(A) + F23_3_t(A))

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------



# ------------------------------------------------------------------------------
# ---- 9. THE SPHERICAL CAPS; S̃ -----------------------------------------------

# The height of the spherical cap (any dimensional spherical cap).
function SphericalCapHeight(H) # H (capital η) is the distance from the origin to the base of the spherical cap

    # Finds the height of the spherical cap as an interval.
    Height = @interval(1 - H)

    # Since the height must be between 0 and 1 inclusive, if the interval for height exceeds those boundaries, then we replace the lower bound for height (which is < 0) by 0 and/or the upper bound for height (which is > 1) by 1.
    height_lo = min(max(0, IntervalArithmetic.inf(Height)), 1)
    height_hi = max(0, min(1, IntervalArithmetic.sup(Height)))
    Height_adjusted = IntervalArithmetic.interval(height_lo, height_hi) # The use of "interval" instead of "@interval" is okay here because the endpoints for Height already contain 1 - H. 

    return Height_adjusted

end

# The surface area of the curved portion of the spherical cap with height η.

# The (2-D) surface area of the curved portion of a spherical cap of B^3.
S123(A) = @interval((2 * π)) * SphericalCapHeight(H123(A))

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------



# ------------------------------------------------------------------------------
# ---- 10. INTEGRALS; ω --------------------------------------------------------

function Integralequation(f, a, b)
    v = quadgk_count(f, a, b)[1]
    e = quadgk_count(f, a, b)[2]
    return I(v - e, v + e)
end 

function IntegralEquation(F_inf, F_sup, A, B)
    V_inf = quadgk_count(F_inf, sup(A), inf(B))[1]
    V_sup = quadgk_count(F_sup, inf(A), sup(B))[1]
    E_inf = quadgk_count(F_inf, sup(A), inf(B))[2]
    E_sup = quadgk_count(F_sup, inf(A), sup(B))[2]
    return I(V_inf - E_inf, V_sup + E_sup)
end 

# Some single letters (e.g. "a") give errors when assigned as names of functions, so I use the underscore to get around this problem.

_k(m, r)    = floor((1 + r^2) / (1 - r^2) + m) + 1 # The function k [page 210]
_a(r)       = √(1 - r^2) # This function is a(r) in Betke and Henk (1998)
_b(k, m)    = √((k + 1 - m) / (k + 3 - m))
μ0(k, m, r) = r / _b(k, m) # The function μ0 

integrand_m(d, l, k, m, r, μ) = (√((_a(r))^2 + μ^2 * (_b(k, m))^2))^(-(d - l + m)) * μ^(d - l + m - (k + 2)) * (1 - μ)^k 
m0(d, l, k, m, r) = Integralequation(μ -> integrand_m(d, l, k, m, r, μ), 0, μ0(k, m, r))
m1(d, l, k, m, r) = Integralequation(μ -> integrand_m(d, l, k, m, r, μ), μ0(k, m, r), 1)

# The interval versions of the above. 
# "MU" is capital μ. 
_K(m, r)     = @interval(floor((1 + r^2) / (1 - r^2) + m) + 1)
_A(r)        = @interval(√(1 - r^2))
_B(k, m)     = @interval(√((k + 1 - m) / (k + 3 - m)))
MU0(k, m, r) = @interval(r / _B(k, m))

Integrand_M(d, l, k, m, r, μ) = @interval( (√((_A(r))^2 + μ^2 * (_B(k, m))^2))^(-(d - l + m)) * μ^(d - l + m - (k + 2)) * (1 - μ)^k )
M0(d, l, k, m, r) = IntegralEquation(μ -> inf(Integrand_M(d, l, k, m, r, μ)), μ -> sup(Integrand_M(d, l, k, m, r, μ)), @interval(0), MU0(k, m, r))
M1(d, l, k, m, r) = IntegralEquation(μ -> inf(Integrand_M(d, l, k, m, r, μ)), μ -> sup(Integrand_M(d, l, k, m, r, μ)), MU0(k, m, r), @interval(1))

# The upper bound for q: q ≤ M0 / M1 
q_UB(d, l, m, r) = m0(d, l, _k(m, r), m, r)[1] / m1(d, l, _k(m, r), m, r)[1]
Q_UB(d, l, m, r) = M0(d, l, _K(m, r), m, r) / M1(d, l, _K(m, r), m, r)

# The integrand in the definition of the ω functions.
#     l = the dimension of the polytope P.
#     m = the dimension of the normal cone.
#     p = power.
# The integrand is bounded below by 1 / (1 + M0 / M1).
# Also, 0^0 and 0.0^0.0 both equal 1 in Julia, so we can write "r^p" without other conditions.
Integrand_Ω(d, l, m, r, p) = inf(@interval(r^p / (1 + Q_UB(d, l, m, r))))
Ω(d, l, m, p, lb, ub) = IntegralEquation(r -> inf(Integrand_Ω(d, l, m, r, p)), r -> sup(Integrand_Ω(d, l, m, r, p)), lb, ub)

u(d) = min(0.95, 1 - 9 / d^(3/2))
U(d) = min(inf(@interval(0.95)), inf(@interval(1 - 9 / d^(3/2))))

Ω1_2(d) = Ω(d, 2, 1, 0, 0, U(d))

Ω1_3(d)         = Ω(d, 3, 2, 0, 0, U(d))
Ω2_3_min(α★, d) = Ω(d, 3, 1, 1, 0, min(inf(H_min_2(α★)), U(d)))
Ω3_3(d)         = Ω(d, 3, 1, 1, 0, U(d))
Ω3_3_min(α★, d) = Ω(d, 3, 0, 2, 0, min(inf(H123_min(α★)), U(d)))

Ω̂3_3(α★, d)     = Ω(d, 3, 0, 2, min(inf(H123_min(α★)), inf(U(d))), U(d))

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------



# ------------------------------------------------------------------------------
# ---- 11. SETUP FOR THE VOLUME LOWER BOUNDS -----------------------------------

# In this section, α can be α21, α31, or α41. 

# The g functions. 

μ(α, ζ) = (√(4 - ζ^2) - 2 * sin(α)) / (2 + ζ - 2 * sin(α))

integrand_g0(m, r,    ζ, d) = r^(m - 2) * (r * ζ / √(4 - ζ^2) + 2 / √(4 - ζ^2))^(d - (m - 1))
integrand_g1(m, r, α, ζ, d) = r^(m - 2) * (r * (sin(α) - 1) / sin(α) * √((2 + ζ) / (2 - ζ)) + 1 / sin(α))^(d - (m - 1))

g0(m, α, ζ, d) = Integralequation(r -> integrand_g0(m, r, ζ, d), 0, μ(α, ζ))
g1(m, α, ζ, d) = Integralequation(r -> integrand_g0(m, r, ζ, d), μ(α, ζ), √((2 - ζ) / (2 + ζ)))
 g(m, α, ζ, d) = g0(m, α, ζ, d)[1] + g1(m, α, ζ, d)[1]

function g_min(m, α, d) 
    
    ζ_max = min(2 * sin(α), 2 * cos(α))
    
    return min(g(m, α, 0, d), g(m, α, ζ_max, d)) 
    # If we check ζ = 0 and ζ = ζ_max = min(2 * sin(α), 2 * cos(α)) then it should be fine. 

end 

MU(α, ζ) = intersect_interval(@interval((√(4 - ζ^2) - 2 * sin(α)) / (2 + ζ - 2 * sin(α))), I(0, 1))

Integrand_G0(m, r,    ζ, d) = @interval( r^(m - 2) * (r * ζ / √(4 - ζ^2) + 2 / √(4 - ζ^2))^(d - (m - 1)) )
Integrand_G1(m, r, α, ζ, d) = @interval( r^(m - 2) * (r * (sin(α) - 1) / sin(α) * √((2 + ζ) / (2 - ζ)) + 1 / sin(α))^(d - (m - 1)) )

G0(m, α, ζ, d) = IntegralEquation(r -> inf(Integrand_G0(m, r, ζ, d)), r -> sup(Integrand_G0(m, r, ζ, d)), @interval(0), MU(α, ζ))
G1(m, α, ζ, d) = IntegralEquation(r -> inf(Integrand_G0(m, r, ζ, d)), r -> sup(Integrand_G0(m, r, ζ, d)), MU(α, ζ), @interval(√((2 - ζ) / (2 + ζ))))
 G(m, α, ζ, d) = G0(m, α, ζ, d)[1] + G1(m, α, ζ, d)[1]

function G_min(m, α, d) 
    
    ζ_max = min(sup(2 * sin(@interval(α))), sup(2 * cos(@interval(α))))
    
    return min(G(m, α, 0, d), G(m, α, ζ_max, d)) 
    # If we check ζ = 0 and ζ = ζ_max = min(2 * sin(α), 2 * cos(α)) then it should be fine. 

end 

# The p functions. 

integrand_p(α, d, r, m) = r^(m - 2) * (-r * cos(α) / sin(α) + 1 / sin(α))^(d - (m - 1)) # FFS I wrote (d - m - 1) instead of (d - (m - 1)) and kept getting wrong answers 
Integrand_P(α, d, r, m) = @interval(r^(m - 2) * (-r * cos(α) / sin(α) + 1 / sin(α))^(d - (m - 1)))

p(m, α, d) = Integralequation(r -> integrand_p(α, d, r, m), 0, (1 - sin(α)) / cos(α))
P(m, α, d) = IntegralEquation(r -> inf(Integrand_P(α, d, r, m)), r -> sup(Integrand_P(α, d, r, m)), @interval(0), @interval((1 - sin(α)) / cos(α)))

# The q functions (p1, p2, etc. in BH1998) that combine g and p. 

function Q(m, α, d) 
    if α < π / 4 
        return G_min(m, α, d)
    else # π / 4 ≤ α ≤ α★[m] 
        return min(G_min(m, α, d), P(m, α, d))
    end 
end 

function Q̂(m, α, d) 
    if α < π / 4 
        return G_min(m, α, d)
    else # π / 4 ≤ α ≤ α★[m] 
        return min(G_min(m, α, d), @interval(2 * P(m, α, d)))
    end 
end 

# The ν function. 
ν(θ) = (π - θ) / 2 - (acos(2 * sin(θ / 2)) - 2 * sin(θ / 2) * √(1 - (2 * sin(θ / 2))^2))
N(Θ) = @interval((π - Θ) / 2 - (acos(2 * sin(Θ / 2)) - 2 * sin(Θ / 2) * √(1 - (2 * sin(Θ / 2))^2)))

function Ω20(α★, d)

    Ω2_part = Ω2_3_min(α★, d)
    K2_part = K(d - 2) / K(d - 1) 
    K0_part = K(d) / K(d - 1)
    Total = @interval(1 / 2 * Ω2_part) * K(d - 2) / K(d - 1) - @interval(1 / (4 * π)) * K(d) / K(d - 1)

    println("    Ω2_3_min(α★, d) = $(inf(Ω2_part))")
    println("K(d - 2) / K(d - 1) = $(inf(K2_part))")
    println("K(d)     / K(d - 1) = $(inf(K0_part))")
    println("")
    println("         Ω20(α★, d) = $Total")

end

function Ω31(α★, d) 

    Ω3_part = Ω3_3_min(α★, d) 
    K3_part = K(d - 3) / K(d - 1) 
    Ω1_part = Ω1_3(d)
    Total = -Ω3_part * K3_part +  @interval(1 / (2 * π)) * Ω1_part

    println("    Ω3_3_min(α★, d) = $(inf(Ω3_part))")
    println("K(d - 3) / K(d - 1) = $(inf(K3_part))")
    println("")
    println("        Ω31(α★, d)  = $Total")

end 

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------



# ------------------------------------------------------------------------------
# ---- 12. THE VOLUME LOWER BOUNDS FOR ALL CASES -------------------------------

#= Key: For the lower bound expressions below, the following suffixes are used. 

"P" = Positive (or zero) y1 ⋅ y2 
"N" = Negative y1 ⋅ y2 

"T" = Tiny angle 
"S" = Small angle 
"M" = Medium angle 
"L" = Large angle 

BHW1994: 
1:    "e"           0 ≤ Â21 ≤ π / 4    0 ≤ y1 ⋅ y2 
2:    "T"           0 ≤ Â21 ≤ π / 4        y1 ⋅ y2 < 0 
3(a): "M"       π / 4 ≤ Â21 ≤ A2⋆ 
3(b): "L"         A2⋆ ≤ Â21 

BH1998: 
1:    "e"           0 ≤ Â21 < π / 3    0 ≤ y1 ⋅ y2 
2(a): "T"           0 ≤ Â21 ≤ π / 4        y1 ⋅ y2 < 0 
2(b): "S"       π / 4 ≤ Â21 < π / 3        y1 ⋅ y2 < 0 
3(a): "MS"      π / 3 ≤ Â21 ≤ A2⋆ 
3(b): "LS"        A2⋆ ≤ Â21                  A31 ≤ A3⋆ 
4:    "LL"     B(A3⋆) ≤ Â21            A3⋆ ≤ A31 
4':   "LL'"    B(A3⋆) ≤ Â21            A3⋆ ≤ A31         (also includes D0) 
=#



# BHW1994: 

# Case 1: Endpoint. 
#         1       
# vol D ≥ ─ * κ(d)
#         2       
#             ⎡ 3 - 2 * √2   κ(d - 2)   κ(d - 1)     1   ⎤ 
#           + ⎢ ────────── * ──────── - ──────── - ───── ⎥ * α̂[2,1] * κ(d)
#             ⎣     2          κ(d)       κ(d)     2 * π ⎦ 
# 
#           + κ(d - 1) 
# 
#         1 
#       = ─ * κ(d) + λ̂_endpoint_BHW1994(d) * α̂[2,1] * κ(d) + κ(d - 1) 
#         2 
#         1 
#       > ─ * κ(d) + κ(d - 1). 
#         2 
Λ̂_endpoint_BHW1994(d) = inf(@interval((3 - 2 * √2) / 2)) * κ_big(d - 2) / κ_big(d) - κ_big(d - 1) / κ_big(d) - inf(@interval(1 / (2 * π))) 

# Case 2: Acute A21. 
# vol D ≥ 2 * κ(d - 1) 
#             ⎡ 3 - 2 * √2   κ(d - 2)   π     ⎞ ⎤ 
#           + ⎢ ────────── * ──────── + ─ - 2 ⎟ ⎥ * α̂[2,1] * κ(d - 1) 
#             ⎣     2        κ(d - 1)   8     ⎠ ⎦ 
#       = 2 * κ(d - 1) + λ̂_A_BHW1994(d) * α̂[2,1] * κ(d) + κ(d - 1) 
#       > 2 * κ(d - 1). 
Λ̂_A_BHW1994(d) = inf( @interval((3 - 2 * √2) / 2) * κ_big(d - 2) / κ_big(d - 1) - @interval(2 + π / 8) )

# Case 3 (a): Medium A21. 
# vol D ≥ 2 * κ(d - 1) 
#             ⎡ α̂[2,1]   1 - sin(α⋆[2])   κ(d - 2)   1 - sin((α⋆[2])     ⎤ 
#           + ⎢ ────── * ────────────── * ──────── + ─────────────── - 2 ⎥ 
#             ⎣   2      1 + sin(α⋆[2])   κ(d - 1)     cos(α⋆[2])        ⎦ 
#       = 2 * κ(d - 1) + λ_M_BHW1994(d) * κ(d - 1) 
#       > 2 * κ(d - 1). 
Λ_M_BHW1994(α★, d) = inf( @interval(α★[2] / 2) * @interval((1 - sin(α★[2])) / (1 + sin(α★[2]))) * κ_big(d - 2) / κ_big(d - 1) + @interval((1 - sin(α★[2])) / cos(α★[2]) - 2) )

# Case 3 (b): Large A21. 
# vol D ≥ 2 * κ(d - 1) 
#          ⎡ α̂[2,1]     ε               1               κ(d - 2) 
#        + ⎢ ────── * ───── * ─────────────────────── * ──────── 
#          ⎣   2      ε + 2   1 + q_BHW1994(1 + ε, d)   κ(d - 1) 
#                                                         1 - sin((α⋆[2])     ⎤ 
#                                                       + ─────────────── - 2 ⎥ 
#                                                           cos(α⋆[2])        ⎦ 
#       = 2 * κ(d - 1) + λ̂_L_BHW1994(d) * κ(d - 1) 
#       > 2 * κ(d - 1). 
integrand_m_BHW1994(x, d) = (1 - BigFloat(x)^2)^((d - 5) / 2)
m0_BHW1994(ε, d) = quadgk_count(x -> integrand_m_BHW1994(x, d), √3 / 2, 1 / ε)[1]
m1_BHW1994(ε, d) = quadgk_count(x -> integrand_m_BHW1994(x, d), 1 / ε, 1)[1]
q_BHW1994(ε, d) = BigFloat(m1_BHW1994(ε, d) / m0_BHW1994(ε, d))
Λ_L_BHW1994(α★, ε, d) = inf(@interval(α★[2] / 2 * ε / (ε + 2) * 1)) / (1 + q_BHW1994(1 + ε, d)) * Float64(κ_big(d - 2) / κ_big(d - 1)) + inf(@interval((1 - sin(α★[2])) / cos(α★[2]) - 2))



# BH1998 and this thesis: 

# Case 1: Endpoint. 
#         1       
# vol D ≥ ─ * κ(d)
#         2       
#             ⎛                  κ(d - 2)   1  ⎞   Â[2,1]
#           + ⎜ q(3, α⋆[3], d) * ──────── - ── ⎟ * ────── * κ(d) + κ(d - 1)
#             ⎝                    κ(d)     2π ⎠     2
#         1                              α̂[2,1]
#       = ─ * κ(d) + λ_endpoint(A⋆, d) * ────── * κ(d) + κ(d - 1).
#         2                                2
Λ_endpoint(α★, d) = inf( Q(3, α★[3], d) * K(d - 2) / K(d) - @interval(1 / (2 * π)) )

# Case 2 (a): Acute A21. 
# vol D ≥ [vol D̄1(P2) + vol D̂1(P2)] + vol D2(P2) 
#       ≥ 2 * κ(d - 1)
#               ⎧    ⎡  ⎛ π ⎞   κ(d - 2)                   ⎤     ⎫ 
#          + min⎨ 0, ⎢ ν⎜ ─ ⎟ * ──────── +  q(2, α⋆[2], d) ⎥ - 2 ⎬ * κ(d - 1) 
#               ⎩    ⎣  ⎝ 8 ⎠   κ(d - 1)                   ⎦     ⎭ 
#       = 2 * κ(d - 1) + min{0, λ_A(d) - 2} * κ(d - 1).  
Λ_A(α★, d) = inf( N(π / 8) * K(d - 2) / K(d - 1) + Q(2, α★[2], d) - @interval(2) )
# Note: Should be vol D ≥ 2 * κ(d - 1) + min{0, λ_A(d)} * κ(d - 1) = 0 but we include this term for extra info. 

# Case 2 (b): Small A21. 
# vol D ≥ vol D0(P2) + vol D̂1(P2) + vol D2(P2) 
#       ≥ 2 * κ(d - 1)
#           ⎛ ⎡      ⎛ α[2,1] ⎞    ⎛   π   ⎞   κ(d - 2)                  ⎤     ⎞
#         + ⎜ ⎢ 2 * ν⎜ ────── ⎟ * q⎜3, ─, d⎟ * ──────── + q(2, α⋆[2], d) ⎥ - 2 ⎟
#           ⎝ ⎣      ⎝   2    ⎠    ⎝   3   ⎠   κ(d - 1)                  ⎦     ⎠
#                                                                    * κ(d - 1) 
#      ≥ 2 * κ(d - 1) + (λ_S(α, d) - 2) * κ(d - 1). 
Λ_S(A★, d, ϕ = π / 3) = inf( @interval(2) * N(π / 8) * Q(3, ϕ, d) * K(d - 2) / K(d - 1) + Q(2, A★[2], d) - @interval(2) )

# Case 3 (a): Medium A21 and small A31. 
# vol D ≥ [vol D0(P2) + vol D2(P2)] + vol D1(P2) 
#       ≥ 2 * κ(d - 1) 
#            ⎛ ⎡                  π                        κ(d - 2) ⎤     ⎞ 
#          + ⎜ ⎢ q(2, α⋆[2], d) + ─ * 2 * q̂(3, α⋆[3], d) * ──────── ⎥ - 2 ⎟ 
#            ⎝ ⎣                  3                        κ(d - 1) ⎦     ⎠ 
#                                                                    * κ(d - 1) 
#       ≥ 2 * κ(d - 1) * (λ_MS(α⋆, d) - 2) * κ(d - 1). 
Λ_MS(A★, d, ϕ = π / 3) = inf( Q(2, A★[2], d) + @interval(ϕ) * Q̂(3, A★[3], d) * K(d - 2) / K(d - 1) - @interval(2) )

# Case 3 (b): Large A21 and small A31. 
# vol D ≥ [vol D0(P2) + vol D2(P2)] + vol D1(P2) 
#       ≥ 2 * κ(d - 1) + (λ02_2(α⋆, d) + λ1_2(α⋆, d)) * κ(d - 1) 
#       ≥ 2 * κ(d - 1)
#            ⎛ ⎡           α⋆[2]                        κ(d - 2) ⎤     ⎞ 
#          + ⎜ ⎢ λ1_2(d) + ───── * 2 * q̂(3, α⋆[3], d) * ──────── ⎥ - 2 ⎟ 
#            ⎝ ⎣             2                          κ(d - 1) ⎦     ⎠ 
#                                                                    * κ(d - 1) 
#       ≥ 2 * κ(d - 1) + (λ_LS(α⋆, d) - 2) * κ(d - 1). 
function Λ_LS(α★, d)

    # Components of the volume lower bound. 
    Λ02_2(α★, d) = @interval(α★[2] / 2 * 2) * Q̂(3, α★[3], d) 
         Λ1_2(d) = Ω1_2(d)

    volumeLowerBound = Λ1_2(d) * K(d - 1) + Λ02_2(α★, d) * K(d - 2)
    
    return inf( volumeLowerBound / K(d - 1) - @interval(2) )

end 

# Case 4: Large A21 and large A31. 

Λ0_3(A, d)        = G123_3(A) / @interval(4 * π) * K(d) / K(d - 1)
Λ1_3(A, d)        = (G23_3_t(A) + G13_3_t(A) + G12_3_t(A)) / @interval(2 * π) * Ω1_3(d)
Λ2_3(A, α★, d)    = (F12_3_t(A) + F13_3_t(A) + F23_3_t(A)) / @interval(2) * Ω2_3_min(α★, d) * K(d - 2) / K(d - 1)
Λ3_3(A, α★, d)    = F123_3(A) * Ω3_3_min(α★, d) * K(d - 3) / K(d - 1)
Λ123_3(A, α★, d)  = [Λ1_3(A, d) Λ2_3(A, α★, d) Λ3_3(A, α★, d)] 
Λ0123_3(A, α★, d) = [Λ0_3(A, d) Λ1_3(A, d) Λ2_3(A, α★, d) Λ3_3(A, α★, d)] 

A_min(α★) = ◺([B(@interval(α★[3]))   @interval(α★[3])   B(@interval(α★[3])) / 2])

# vol D ≥ [vol D1(P3) + vol D2(P3)] + vol D3(P3) 
#       ≥ 2 * κ(d - 1) + λ123_3(α, d) * κ(d - 1) 
#       ≥ 2 * κ(d - 1) 
#            ⎛       ⎛ ⎛                  β(α⋆[3]) ⎞    ⎞     ⎞ 
#          + ⎜ λ123_3⎜ ⎜ β(α⋆[3]), α⋆[3], ──────── ⎟, d ⎟ - 2 ⎟ * κ(d - 1)
#            ⎝       ⎝ ⎝                     2     ⎠    ⎠     ⎠
#       ≥ 2 * κ(d - 1) + (λ_LL(α⋆, d) - 2) * κ(d - 1). 

Λ_LL(α★, d) = inf.(sum(Λ123_3(A_min(α★), α★, d)) - @interval(2))

# Case 4': Large A21 and large A31, and including D0(P3)
# vol D ≥ vol D0(P3) + [vol D1(P3) + vol D2(P3)] + vol D3(P3) 
#       ≥ 2 * κ(d - 1) + λ0123_3(α, d) * κ(d - 1) 
#       ≥ 2 * κ(d - 1) 
#            ⎛        ⎛ ⎛                  β(α⋆[3]) ⎞    ⎞     ⎞ 
#          + ⎜ λ0123_3⎜ ⎜ β(α⋆[3]), α⋆[3], ──────── ⎟, d ⎟ - 2 ⎟ * κ(d - 1)
#            ⎝        ⎝ ⎝                     2     ⎠    ⎠     ⎠
#       ≥ 2 * κ(d - 1) + (λ_LL0(α⋆, d) - 2) * κ(d - 1). 

Λ_LL0(α★, d) = inf.(sum(Λ0123_3(A_min(α★), α★, d)) - @interval(2))

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------



# ------------------------------------------------------------------------------
# ---- # 13. Different numbers of endpoints ------------------------------------

Λ_min_BHW1994(α★, ε, d) = min(Λ_M_BHW1994(α★, d), Λ_L_BHW1994(α★, ε, d))

B_BHW1994(α★, ε, d) = (1 - 1 / 2 * κ_big(d) / κ_big(d - 1)) / (-Λ_min_BHW1994(α★, ε, d)) 

# The Sausage Conjecture holds for these dimensions and numbers of spheres: 
endpoints_2_BHW1994(α★, ε, d) = 2
endpoints_1_BHW1994(α★, ε, d) = Float64(1 + inf(B_BHW1994(α★, ε, d)))
endpoints_0_BHW1994(α★, ε, d) = Float64(2 * inf(B_BHW1994(α★, ε, d)))

Λ_min_BH1998(α★, d) = min(Λ_S(α★, d, π / 3), Λ_MS(α★, d, π / 3), Λ_LS(α★, d), Λ_LL(α★, d))

B_BH1998(α★, d) = (1 - 1 / 2 * K(d) / K(d - 1)) / (-Λ_min(α★, d)) 

# The Sausage Conjecture holds for these dimensions and numbers of spheres: 
endpoints_2_BH1998(α★, d) = 2
endpoints_1_BH1998(α★, d) = 1 + inf(B_BH1998(α★, d))
endpoints_0_BH1998(α★, d) = 2 * inf(B_BH1998(α★, d))

Λ_min(α★, d) = min(Λ_S(α★, d, π / 3), Λ_MS(α★, d, π / 3), Λ_LS(α★, d), Λ_LL0(α★, d))

B(α★, d) = (1 - 1 / 2 * K(d) / K(d - 1)) / (-Λ_min(α★, d)) 

# The Sausage Conjecture holds for these dimensions and numbers of spheres: 
endpoints_2(α★, d) = 2
endpoints_1(α★, d) = 1 + inf(B(α★, d))
endpoints_0(α★, d) = 2 * inf(B(α★, d))

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------



# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

function volumeLowerBoundList(paper, d, ε, α★, ϕ = π / 3)
    
    println("")
    println("Excess volumes:")
    println("")

    if paper == "BHW1994" 

        caseVector = [Λ̂_endpoint_BHW1994(d) Λ̂_A_BHW1994(d) Λ_M_BHW1994(α★, d) Λ_L_BHW1994(α★, ε, d)]
        
        println("                      α⋆[2]")
        println("    Case 1:     Λ̂_e           $(Float64(caseVector[1]))") 
        println("    Case 2:     Λ̂_A     ↘     $(Float64(caseVector[2]))") 
        println("    Case 3 (a): Λ_M     ↘     $(Float64(caseVector[3]))") 
        println("    Case 3 (b): Λ_L     ↗     $(Float64(caseVector[4]))") 
        println("")

        min_Λ = Float64(min(caseVector[3], caseVector[4]))
        println("    Minimum:                  $min_Λ")
        println("")
        
        println("Endpoints:")
        println("")
        if min_Λ < 0 
            println("    2 endpoints: n ≤ $(endpoints_2_BHW1994(α★, ε, d))")
            println("    1 endpoint : n ≤ $(endpoints_1_BHW1994(α★, ε, d))")
            println("    0 endpoints: n ≤ $(endpoints_0_BHW1994(α★, ε, d))")
        else #min_Λ ≥ 0
            println("    2 endpoints: n < ∞")
            println("    1 endpoint : n < ∞")
            println("    0 endpoints: n < ∞")
        end 

    elseif paper == "H1994" 
        return 0 
    elseif paper == "BH1998" 

        caseVector = [Λ_endpoint(α★, d) Λ_A(α★, d) Λ_S(α★, d) Λ_MS(α★, d) Λ_LS(α★, d) Λ_LL(α★, d)]

        println("                      α⋆[2] α⋆[3]")
        println("    Case 1:     Λ̂_e                $(caseVector[1])    (excess)")
        println("    Case 2 (a): Λ_A     ↘     -    $(caseVector[2])") 
        println("    Case 2 (b): Λ_S     ↘     -    $(caseVector[3])") 
        println("    Case 3 (a): Λ_MS    ↘     ↘    $(caseVector[4])") 
        println("    Case 3 (b): Λ_LS    ↗     ↘    $(caseVector[5])") 
        println("    Case 4:     Λ_LL    -     ↗    $(caseVector[6])") 
        println("                   vol(D^1(P^3)) ≥ $(inf(Λ1_3(A_min(α★), d)))") 
        println("                   vol(D^2(P^3)) ≥ $(inf(Λ2_3(A_min(α★), α★, d)))") 
        println("                   vol(D^3(P^3)) ≥ $(inf(Λ3_3(A_min(α★), α★, d)))") 
        println("")
        
        min_Λ = min(caseVector[2], caseVector[3], caseVector[4], caseVector[5], caseVector[6])
        println("    Minimum:                       $min_Λ")
        println("")
        
        println("Endpoints:")
        println("")
        if min_Λ < 0 
            println("    2 endpoints: n ≤ $(endpoints_2_BH1998(α★, d))")
            println("    1 endpoint : n ≤ $(endpoints_1_BH1998(α★, d))")
            println("    0 endpoints: n ≤ $(endpoints_0_BH1998(α★, d))")
        else #min_Λ ≥ 0
            println("    2 endpoints: n < ∞")
            println("    1 endpoint : n < ∞")
            println("    0 endpoints: n < ∞")
        end 

    else #if paper == "Thesis"

        caseVector = [Λ_endpoint(α★, d) Λ_A(α★, d) Λ_S(α★, d, ϕ) Λ_MS(α★, d, ϕ) Λ_LS(α★, d) Λ_LL0(α★, d)]

        println("                      α⋆[2] α⋆[3]")
        println("    Case 1:     Λ̂_e                $(caseVector[1])    (excess)")
        println("    Case 2 (a): Λ_A     ↘     -    $(caseVector[2])") 
        println("    Case 2 (b): Λ_S     ↘     -    $(caseVector[3])") 
        println("    Case 3 (a): Λ_MS    ↘     ↘    $(caseVector[4])") 
        println("    Case 3 (b): Λ_LS    ↗     ↘    $(caseVector[5])") 
        println("    Case 4:     Λ_LL0   -     ↗    $(caseVector[6])") 
        println("")
        println("                   vol(D^0(P^3)) ≥ $(inf(Λ0_3(A_min(α★), d)))") 
        println("                   vol(D^1(P^3)) ≥ $(inf(Λ1_3(A_min(α★), d)))") 
        println("                   vol(D^2(P^3)) ≥ $(inf(Λ2_3(A_min(α★), α★, d)))") 
        println("                   vol(D^3(P^3)) ≥ $(inf(Λ3_3(A_min(α★), α★, d)))") 
        println("")
        
        min_Λ = min(caseVector[2], caseVector[3], caseVector[4], caseVector[5], caseVector[6])
        println("    Minimum:                       $min_Λ")
        println("")
        
        println("Endpoints:")
        println("")
        if min_Λ < 0 
            println("    2 endpoints: n ≤ $(endpoints_2(α★, d))")
            println("    1 endpoint : n ≤ $(endpoints_1(α★, d))")
            println("    0 endpoints: n ≤ $(endpoints_0(α★, d))")
        else #min_Λ ≥ 0
            println("    2 endpoints: n < ∞")
            println("    1 endpoint : n < ∞")
            println("    0 endpoints: n < ∞")
        end 

    end 

end 

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------