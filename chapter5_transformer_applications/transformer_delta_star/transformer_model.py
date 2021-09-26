import math

dt = 1.0e-6

# Phase voltage rating of windings
# Divide by sqrt(3) for star connection
Vrated_primary = 11000.0
Vrated_secondary = 415.0 / math.sqrt(3.0)

# Phase a - Windings 1 and 2
# Phase b - Windings 3 and 4
# Phase c - Windings 5 and 6

# Per-phase power rating
VArated = 50000.0 / 3    # Volt-Amperes
# Primary windings
V1rated = Vrated_primary
V3rated = Vrated_primary
V5rated = Vrated_primary
# Secondary windings
V2rated = Vrated_secondary
V4rated = Vrated_secondary
V6rated = Vrated_secondary
# Rated frequency
frated = 50.0       # Hz
omega_rated = 2*math.pi*frated

# Computations performed at the start of the simulation
if init_simulation > 0:

    #Current rating
    I1rated = VArated / V1rated
    I2rated = VArated / V2rated
    I3rated = VArated / V3rated
    I4rated = VArated / V4rated
    I5rated = VArated / V5rated
    I6rated = VArated / V6rated

    # Rated (base) impedance
    Z1rated = V1rated / I1rated
    Z2rated = V2rated / I2rated
    Z3rated = V3rated / I3rated
    Z4rated = V4rated / I4rated
    Z5rated = V5rated / I5rated
    Z6rated = V6rated / I6rated

    # Self inductance
    L1 = 50.0 * Z1rated / omega_rated
    L2 = 50.0 * Z2rated / omega_rated
    L3 = 50.0 * Z3rated / omega_rated
    L4 = 50.0 * Z4rated / omega_rated
    L5 = 50.0 * Z5rated / omega_rated
    L6 = 50.0 * Z6rated / omega_rated

    # Leakage inductance
    Ll1 = 0.02 * Z1rated / omega_rated
    Ll2 = 0.02 * Z2rated / omega_rated
    Ll3 = 0.02 * Z3rated / omega_rated
    Ll4 = 0.02 * Z4rated / omega_rated
    Ll5 = 0.02 * Z5rated / omega_rated
    Ll6 = 0.02 * Z6rated / omega_rated

    # Magnetizing inductance
    Lm1 = L1 - Ll1
    Lm2 = L2 - Ll2
    Lm3 = L3 - Ll3
    Lm4 = L4 - Ll4
    Lm5 = L5 - Ll5
    Lm6 = L6 - Ll6

    # Mutual inductance within windings of a phase
    M12 = math.sqrt(Lm1 * Lm2)
    M34 = math.sqrt(Lm3 * Lm4)
    M56 = math.sqrt(Lm5 * Lm6)

    # Mutual inductance between phases
    # Half of the flux generated in one limb links with a winding in another limb
    k_factor = 0.5
    M13 = k_factor * math.sqrt(Lm1 * Lm3)
    M14 = k_factor * math.sqrt(Lm1 * Lm4)
    M15 = k_factor * math.sqrt(Lm1 * Lm5)
    M16 = k_factor * math.sqrt(Lm1 * Lm6)
    M23 = k_factor * math.sqrt(Lm2 * Lm3)
    M24 = k_factor * math.sqrt(Lm2 * Lm4)
    M25 = k_factor * math.sqrt(Lm2 * Lm5)
    M26 = k_factor * math.sqrt(Lm2 * Lm6)
    M35 = k_factor * math.sqrt(Lm3 * Lm5)
    M36 = k_factor * math.sqrt(Lm3 * Lm6)
    M45 = k_factor * math.sqrt(Lm4 * Lm5)
    M46 = k_factor * math.sqrt(Lm4 * Lm6)

    # Winding resistance
    r1 = 0.01 * Z1rated
    r2 = 0.01 * Z2rated
    r3 = 0.01 * Z3rated
    r4 = 0.01 * Z4rated
    r5 = 0.01 * Z5rated
    r6 = 0.01 * Z6rated

    # Core loss resistance
    Rc1 = V1rated * V1rated / (0.01 *VArated)
    Rc2 = V2rated * V2rated / (0.01 *VArated)
    Rc3 = V3rated * V3rated / (0.01 *VArated)
    Rc4 = V4rated * V4rated / (0.01 *VArated)
    Rc5 = V5rated * V5rated / (0.01 *VArated)
    Rc6 = V6rated * V6rated / (0.01 *VArated)

    # Resistance of VariableResistors in series with
    # ControlledVoltageSources
    res_output1 = 100.0 * Z1rated
    res_output2 = 100.0 * Z2rated
    res_output3 = 100.0 * Z3rated
    res_output4 = 100.0 * Z4rated
    res_output5 = 100.0 * Z5rated
    res_output6 = 100.0 * Z6rated

    # Initializing winding currents at start of simulation
    winding_currents = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

    L = [
        [L1, M12, M13, M14, M15, M16],
        [M12, L2, M23, M24, M25, M26],
        [M13, M23, L3, M34, M35, M36],
        [M14, M24, M34, L4, M45, M46],
        [M15, M25, M35, M45, L5, M56],
        [M16, M26, M36, M46, M56, L6]
    ]
    R = [
        [r1, 0.0, 0.0, 0.0, 0.0, 0.0],
        [0.0, r2, 0.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, r3, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, r4, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.0, r5, 0.0],
        [0.0, 0.0, 0.0, 0.0, 0.0, r6]
    ]
    B = [
        [1.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [0.0, 1.0, 0.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 1.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 1.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.0, 1.0, 0.0],
        [0.0, 0.0, 0.0, 0.0, 0.0, 1.0]
    ]

    # Triangularization
    for count1 in range(len(L)):
        if not L[count1][count1]:
            for count2 in range(count1+1, len(L)):
                if L[count2][count1]:
                    L[count1], L[count2] = L[count2], L[count1]
                    R[count1], R[count2] = R[count2], R[count1]
                    B[count1], B[count2] = B[count2], B[count1]
                    break

        if L[count1][count1]:
            for count2 in range(count1+1, len(L)):
                comm_factor = L[count2][count1]/L[count1][count1]
                for count3 in range(len(L[count1])):
                    L[count2][count3] -= L[count1][count3]*comm_factor
                    R[count2][count3] -= R[count1][count3]*comm_factor
                    B[count2][count3] -= B[count1][count3]*comm_factor

    # End of initialization by setting flag negative
    init_simulation = -1


if t_clock >= t1:
    # Model

    V = [v1, v2, v3, v4, v5, v6]

    # Numerical integration
    dibydt = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    for count1 in range(len(L)-1, -1, -1):
        # If the diagonal element of L is non-zero
        # The row is a differential equation
        if L[count1][count1]:
            k = [0.0, 0.0, 0.0, 0.0]
            for k_count in range(len(k)):
                k[k_count] = 0.0
                for count2 in range(len(B)):
                    k[k_count] += B[count1][count2]*V[count2]
                for count2 in range(count1+1, len(L)):
                    k[k_count] -= L[count1][count2]*dibydt[count2]
                for count2 in range(len(R)):
                    if k_count==0:
                        if count2>count1:
                            k[k_count] -= R[count1][count2]*(winding_currents[count2] + dibydt[count2]*dt)
                        else:
                            k[k_count] -= R[count1][count2]*winding_currents[count2]
                    if k_count==1:
                        if count2==count1:
                            k[k_count] -= R[count1][count2]*(winding_currents[count2] + k[0]*dt/2.0)
                        elif count2>count1:
                            k[k_count] -= R[count1][count2]*(winding_currents[count2] + dibydt[count2]*dt)
                        else:
                            k[k_count] -= R[count1][count2]*winding_currents[count2]
                    if k_count==2:
                        if count2==count1:
                            k[k_count] -= R[count1][count2]*(winding_currents[count2] + k[1]*dt/2.0)
                        elif count2>count1:
                            k[k_count] -= R[count1][count2]*(winding_currents[count2] + dibydt[count2]*dt)
                        else:
                            k[k_count] -= R[count1][count2]*winding_currents[count2]
                    if k_count==3:
                        if count2==count1:
                            k[k_count] -= R[count1][count2]*(winding_currents[count2] + k[2]*dt)
                        elif count2>count1:
                            k[k_count] -= R[count1][count2]*(winding_currents[count2] + dibydt[count2]*dt)
                        else:
                            k[k_count] -= R[count1][count2]*winding_currents[count2]
                k[k_count] = k[k_count]/L[count1][count1]
            dibydt[count1] = (k[0] + 2.0*k[1] + 2.0*k[2] + k[3])/6.0
            winding_currents[count1] += dibydt[count1]*dt
        # If the diagonal element of L is zero
        # The row is an algebraic equation
        else:
            winding_currents[count1] = 0.0
            for count2 in range(len(B)):
                winding_currents[count1] += B[count1][count2]*V[count2]
            for count2 in range(count1+1, len(L)):
                winding_currents[count1] -= L[count1][count2]*dibydt[count2]
            for count2 in range(count1+1, len(R)):
                if not count2==count1:
                    winding_currents[count1] -= R[count1][count2]*winding_currents[count2]
            winding_currents[count1] = winding_currents[count1] / R[count1][count1]

    vout1 = v1 - winding_currents[0]*res_output1
    vout2 = v2 - winding_currents[1]*res_output2
    vout3 = v3 - winding_currents[2]*res_output3
    vout4 = v4 - winding_currents[3]*res_output4
    vout5 = v5 - winding_currents[4]*res_output5
    vout6 = v6 - winding_currents[5]*res_output6

    transf1_circ_curr = iprim1 + iprim2 + iprim3

    t1 += dt
