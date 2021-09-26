import math

dt = 1.0e-7

# Inductances
L1 = 0.001
L2 = 0.001
# Winding resistances
r1 = 0.1
r2 = 0.1

# Coupling factor and mutual inductance
coupling_factor = 0.99
M = coupling_factor*math.sqrt(L1 * L2)

# Resistance of VariableResistors in series with
# ControlledVoltageSources
res_output1 = 100.0
res_output2 = 100.0

if t_clock >= t1:
    # Model
    # e = e - ir
    e1 = v1 - ind_curr1*r1
    e2 = v2 - ind_curr2*r2

    # psi = int e dt
    psi1 += e1*dt
    psi2 += e2*dt

    # psi1 = L1 i1 + M i2
    # psi2 = M i1 + L2 i2
    ind_curr2 = ( psi2 - (M * psi1 / L1 ) ) / ( L2 - (M*M/L1) )
    ind_curr1 = ( psi1 - (M * ind_curr2) ) / L1

    vout1 = v1 - ind_curr1*res_output1
    vout2 = v2 - ind_curr2*res_output2

    coupcoil_e1 = e1
    coupcoil_e2 = e2

    t1 += dt
