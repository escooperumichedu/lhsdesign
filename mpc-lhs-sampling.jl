using SymPy
using DataFrames
using CSV
using Random

T_f = 359.1 # [=] K
x_Af = 1 # unitless
x_Bf = 0 # unitless

# Physical properties
m = 2.79 # [=] mol/kg
rho = 1000 # [=] kg/m3
cp = 4.2 # [=] kJ/(kg*K)
R = 8.314e-3 # [=] kJ/(mol*K)

# Reaction rate constants
k_1 = 9.97e6 # [=] 1/hr
k_2 = 9e6 # [=] 1/hr

# Activation energies and enthalpies of reaction
E_1 = 50 # [=] kJ/mol
E_2 = 60 # [=] kJ/mol
H_1 = -60 # [=] kJ/mol
H_2 = -70 # [=] kJ/mol


# Activity of species
alpha_A = 5 # unitless
alpha_B = 1 # unitless
alpha_C = 0.5 # unitless
epsilon = 0.02


function calculateSS(V1, V2, V3, T1, T2, T3, xA1, xB1, xA2, xB2, xA3, xB3)
    k11 = k_1 * exp(-E_1 / (R * T1))
    k21 = k_2 * exp(-E_2 / (R * T1))
    k12 = k_1 * exp(-E_1 / (R * T2))
    k22 = k_2 * exp(-E_2 / (R * T2))

    # Symbolic variables for flows (input variable)
    @syms Ff1 Ff2 FR F1 F2 F3

    xAR = (alpha_A * xA3) / (alpha_A * xA3 + alpha_B * xB3 + alpha_C * (1 - xB3 - xA3))
    xBR = (alpha_B * xB3) / (alpha_A * xA3 + alpha_B * xB3 + alpha_C * (1 - xB3 - xA3))

    eqn1 = 0 ~ (Ff1 / V1) * (x_Af - xA1) + (FR / V1) * (xAR - xA1) - k11 * xA1 # dxA1/dt
    eqn2 = 0 ~ (F1 / V2) * (xA1 - xA2) + (Ff2 / V2) * (x_Af - xA2) - k12 * xA2 # dxA2/dt
    eqn3 = 0 ~ (F2 / V3) * (xA2 - xA3) - (((0.02 * FR) + FR) / V3) * (xAR - xA3) # dxA3/dt
    eqn4 = 0 ~ Ff1 + FR - F1 # dV1/dt
    eqn5 = 0 ~ Ff2 + F1 - F2 # dV2/dt
    eqn6 = 0 ~ F2 - (0.02 * FR) - FR - F3 # dV3/dt

    # Solve system under constraints
    sol = solve([eqn1, eqn2, eqn3, eqn4, eqn5, eqn6], [F1, F2, F3, Ff1, Ff2, FR])

    F1 = sol[F1]
    F2 = sol[F2]
    F3 = sol[F3]
    FR = sol[FR]
    Ff1 = sol[Ff1]
    Ff2 = sol[Ff2]

    # Symbolic variables for heat rates (input variables)
    @syms Q1 Q2 Q3
    eqn7 = (Ff1 / V1) * (T_f - T1) + (FR / V1) * (T3 - T1) + Q1 / (rho * cp * V1) - (m / cp) * H_1 * k11 * xA1 - (m / cp) * H_2 * k21 * xB1 ~ 0
    eqn8 = (F1 / V2) * (T1 - T2) + (Ff2 / V2) * (T_f - T2) + Q2 / (rho * cp * V2) - (m / cp) * H_1 * k12 * xA2 - (m / cp) * H_2 * k22 * xB2 ~ 0
    eqn9 = (F2 / V3) * (T2 - T3) + Q3 / (rho * cp * V3) ~ 0

    sol = solve([eqn7 eqn8 eqn9], [Q1, Q2, Q3])

    Q1 = sol[Q1]
    Q2 = sol[Q2]
    Q3 = sol[Q3]

    return Float64[F1, F2, F3, Ff1, Ff2, FR, Q1, Q2, Q3]



end

# The origin of which I want to latinhypercube sample perturbations around
V1_sp = 30
V2_sp = 30
V3_sp = 30
T1_sp = 400
T2_sp = 413
T3_sp = 421
xA1_sp = 0.6
xB1_sp = 0.39
xA2_sp = 0.39
xB2_sp = 0.57
xA3_sp = 0.26
xB3_sp = 0.71

x = [V1_sp, V2_sp, V3_sp, T1_sp, T2_sp, T3_sp, xA1_sp, xB1_sp, xA2_sp, xB2_sp, xA3_sp, xB3_sp]
u = (calculateSS(V1_sp, V2_sp, V3_sp, T1_sp, T2_sp, T3_sp, xA1_sp, xB1_sp, xA2_sp, xB2_sp, xA3_sp, xB3_sp))
var = vcat(x, u)'

# Create a dataframe with titles V1 V2 V3 T1 T2 T3 xA1 xB1 xA2 xB2 xA3 xB3 F1 F2 F3 Ff1 Ff2 FR Q1 Q2 Q3
titles = ["V1", "V2", "V3", "T1", "T2", "T3", "xA1", "xB1", "xA2", "xB2", "xA3", "xB3", "F1", "F2", "F3", "Ff1", "Ff2", "FR", "Q1", "Q2", "Q3"]

# Create DataFrame
df = DataFrame(var, Symbol.(titles))
# push!(df, zeros(21))


# Maximum range on vol and temp for changes
V_change = 5
T_change = 5
x_change = 0.10

function lhsdesign(n::Int, p::Int)
    # returns a Latin hypercube sample matrix of size n-by-p. 
    # For each column of X, the n values are randomly distributed with one from each interval (0,1/n), (1/n,2/n), ..., (1 - 1/n,1), and randomly permuted.

    # Empty parameter space
    X = zeros(n, p)

    # n - samples
    # p - variables

    # Window intervals
    a = (0:1:n-1) / n # Left hand side
    b = (1:1:n) / n # Right hand side


    for j in 1:p
        for i in 1:n
            X[i, j] = a[i] + (b[i] - a[i]) * rand()
        end
        X[:, j] .= shuffle!(X[:, j])
    end

    return X
end

# Generate LHS for n samples and p variables
# X = lhsdesign(n, p)

dp = 300 # Number of datapoints we're generating

V1 = zeros(dp)
V2 = zeros(dp)
V3 = zeros(dp)
T1 = zeros(dp)
T2 = zeros(dp)
T3 = zeros(dp)
xA1 = zeros(dp)
xB1 = zeros(dp)
xA2 = zeros(dp)
xB2 = zeros(dp)
xA3 = zeros(dp)
xB3 = zeros(dp)

for k = 1:dp
    println(k)
    X = lhsdesign(1, 9)

    # Print each element
    for i in 1:size(X, 1)
        for j in 1:size(X, 2)
            X[i, j] = rand([-1, 1]) * X[i, j]
        end
    end

    V1[k] = V1_sp + V_change * X[1, 1]
    V2[k] = V2_sp + V_change * X[1, 2]
    V3[k] = V3_sp + V_change * X[1, 3]
    T1[k] = T1_sp + T_change * X[1, 4]
    T2[k] = T2_sp + T_change * X[1, 5]
    T3[k] = T3_sp + T_change * X[1, 6]
    xA1[k] = xA1_sp + x_change * X[1, 7]
    xB1[k] = xB1_sp - x_change * X[1, 7]
    xA2[k] = xA2_sp + x_change * X[1, 8]
    xB2[k] = xB2_sp - x_change * X[1, 8]
    xA3[k] = xA3_sp + x_change * X[1, 9]
    xB3[k] = xB3_sp - x_change * X[1, 9]
    x = [V1[k], V2[k], V3[k], T1[k], T2[k], T3[k], xA1[k], xB1[k], xA2[k], xB2[k], xA3[k], xB3[k]]
    u = (calculateSS(V1[k], V2[k], V3[k], T1[k], T2[k], T3[k], xA1[k], xB1[k], xA2[k], xB2[k], xA3[k], xB3[k]))
    k_vals = vcat(x, u)'
    push!(df, k_vals)
end

# Remove rows that contain any negative elements
filtered_df = filter(row -> all(collect(row) .>= 0), eachrow(df))

using CSV
CSV.write("lhs_sample_data.csv", filtered_df)
