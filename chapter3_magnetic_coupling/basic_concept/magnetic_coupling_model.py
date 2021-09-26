import math

dt = 1.0e-6

# Self-inductances
L1 = 0.3
L2 = 0.3
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

# Initializing inductor currents at start of simulation
if (t_clock <= dt):
    ind_current = [0.0, 0.0]

if t_clock >= t1:
    # Model
    # Matrix representation
    L = [[L1, M], [M, L2]]
    R = [[r1, 0.0], [0.0, r2]]
    V = [v1, v2]

    # Triangularization
    for count1 in range(len(L)):
        if not L[count1][count1]:
            for count2 in range(count1+1, len(L)):
                if L[count2][count1]:
                    L[count1], L[count2] = L[count2], L[count1]
                    R[count1], R[count2] = R[count2], R[count1]
                    V[count1], V[count2] = V[count2], V[count1]
                    break

        if L[count1][count1]:
            for count2 in range(count1+1, len(L)):
                comm_factor = L[count2][count1]/L[count1][count1]
                for count3 in range(len(L[count1])):
                    L[count2][count3] -= L[count1][count3]*comm_factor
                    R[count2][count3] -= R[count1][count3]*comm_factor
                V[count2] -= V[count1]*comm_factor

    # Numerical integration
    dibydt = [0.0, 0.0]
    for count1 in range(len(L)-1, -1, -1):
        # If the diagonal element of L is non-zero
        # The row is a differential equation
        if L[count1][count1]:
            k = [0.0, 0.0, 0.0, 0.0]
            for k_count in range(len(k)):
                k[k_count] = V[count1]
                for count2 in range(count1+1, len(L)):
                    k[k_count] -= L[count1][count2]*dibydt[count2]
                for count2 in range(len(R)):
                    if k_count==0:
                        if count2>count1:
                            k[k_count] -= R[count1][count2]*(ind_current[count2] + dibydt[count2]*dt)
                        else:
                            k[k_count] -= R[count1][count2]*ind_current[count2]
                    if k_count==1:
                        if count2==count1:
                            k[k_count] -= R[count1][count2]*(ind_current[count2] + k[0]*dt/2.0)
                        elif count2>count1:
                            k[k_count] -= R[count1][count2]*(ind_current[count2] + dibydt[count2]*dt)
                        else:
                            k[k_count] -= R[count1][count2]*ind_current[count2]
                    if k_count==2:
                        if count2==count1:
                            k[k_count] -= R[count1][count2]*(ind_current[count2] + k[1]*dt/2.0)
                        elif count2>count1:
                            k[k_count] -= R[count1][count2]*(ind_current[count2] + dibydt[count2]*dt)
                        else:
                            k[k_count] -= R[count1][count2]*ind_current[count2]
                    if k_count==3:
                        if count2==count1:
                            k[k_count] -= R[count1][count2]*(ind_current[count2] + k[2]*dt)
                        elif count2>count1:
                            k[k_count] -= R[count1][count2]*(ind_current[count2] + dibydt[count2]*dt)
                        else:
                            k[k_count] -= R[count1][count2]*ind_current[count2]
                k[k_count] = k[k_count]/L[count1][count1]
            dibydt[count1] = (k[0] + 2.0*k[1] + 2.0*k[2] + k[3])/6.0
            ind_current[count1] += dibydt[count1]*dt
        # If the diagonal element of L is zero
        # The row is an algebraic equation
        else:
            ind_current[count1] = V[count1]
            for count2 in range(count1+1, len(L)):
                ind_current[count1] -= L[count1][count2]*dibydt[count2]
            for count2 in range(count1+1, len(R)):
                if not count2==count1:
                    ind_current[count1] -= R[count1][count2]*ind_current[count2]
            ind_current[count1] = ind_current[count1] / R[count1][count1]

    vout1 = v1 - ind_current[0]*res_output1
    vout2 = v2 - ind_current[1]*res_output2

    e1 = v1 - ind_current[0]*r1
    e2 = v2 - ind_current[1]*r2

    coupcoil_e1 = e1
    coupcoil_e2 = e2

    t1 += dt
