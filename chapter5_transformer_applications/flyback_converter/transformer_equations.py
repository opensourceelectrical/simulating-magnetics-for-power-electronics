import math

dt = 1.0e-8

# Power rating
VArated = 500.0    # Volt-Amperes
# Voltage rating of windings
V1rated = 24.0
V2rated = 15.0
V3rated = 9.0
V4rated = 36.0
# Rated frequency
frated = 10000.0       # Hz
omega_rated = 2*math.pi*frated

# Computations performed at the start of the simulation
if init_simulation > 0:

    #Current rating
    I1rated = VArated / V1rated
    I2rated = VArated / V2rated
    I3rated = VArated / V3rated
    I4rated = VArated / V4rated

    # Rated (base) impedance
    Z1rated = V1rated / I1rated
    Z2rated = V2rated / I2rated
    Z3rated = V3rated / I3rated
    Z4rated = V4rated / I4rated

    # Self inductance
    L1 = 50.0 * Z1rated / omega_rated
    L2 = 50.0 * Z2rated / omega_rated
    L3 = 50.0 * Z3rated / omega_rated
    L4 = 50.0 * Z4rated / omega_rated

    # Leakage inductance
    Ll1 = 0.02 * Z1rated / omega_rated
    Ll2 = 0.02 * Z2rated / omega_rated
    Ll3 = 0.02 * Z3rated / omega_rated
    Ll4 = 0.02 * Z4rated / omega_rated

    # Magnetizing inductance
    Lm1 = L1 - Ll1
    Lm2 = L2 - Ll2
    Lm3 = L3 - Ll3
    Lm4 = L4 - Ll4

    # Mutual inductance
    M12 = math.sqrt(Lm1 * Lm2)
    M13 = math.sqrt(Lm1 * Lm3)
    M14 = math.sqrt(Lm1 * Lm4)
    M23 = math.sqrt(Lm2 * Lm3)
    M24 = math.sqrt(Lm2 * Lm4)
    M34 = math.sqrt(Lm3 * Lm4)

    # Winding resistance
    r1 = 0.01 * Z1rated
    r2 = 0.01 * Z2rated
    r3 = 0.01 * Z3rated
    r4 = 0.01 * Z4rated

    # Core loss resistance
    Rc1 = V1rated * V1rated / (0.01 *VArated)
    Rc2 = V2rated * V2rated / (0.01 *VArated)
    Rc3 = V3rated * V3rated / (0.01 *VArated)
    Rc4 = V4rated * V4rated / (0.01 *VArated)

    # Resistance of VariableResistors in series with
    # ControlledVoltageSources
    srcres_wdg1 = 100.0 * Z1rated
    srcres_wdg2 = 100.0 * Z2rated
    srcres_wdg3 = 100.0 * Z3rated
    srcres_wdg4 = 100.0 * Z4rated

    # Initializing winding currents at start of simulation
    winding_currents = [0.0, 0.0, 0.0, 0.0]

    L = [
        [L1, M12, M13, M14],
        [M12, L2, M23, M24],
        [M13, M23, L3, M34],
        [M14, M24, M34, L4]
    ]
    R = [
        [r1, 0.0, 0.0, 0.0],
        [0.0, r2, 0.0, 0.0],
        [0.0, 0.0, r3, 0.0],
        [0.0, 0.0, 0.0, r4]
    ]
    B = [
        [1.0, 0.0, 0.0, 0.0],
        [0.0, 1.0, 0.0, 0.0],
        [0.0, 0.0, 1.0, 0.0],
        [0.0, 0.0, 0.0, 1.0]
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

    V = [v1meas, v2meas, v3meas, v4meas]

    # Numerical integration
    dibydt = [0.0, 0.0, 0.0, 0.0]
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

    volt_wdg1 = v1meas - winding_currents[0]*srcres_wdg1
    volt_wdg2 = v2meas - winding_currents[1]*srcres_wdg2
    volt_wdg3 = v3meas - winding_currents[2]*srcres_wdg3
    volt_wdg4 = v4meas - winding_currents[3]*srcres_wdg4

    transf1_mag_curr = winding_currents[0] + winding_currents[1]

    transf1_fluxlinkage1 = L[0][0]*winding_currents[0] + L[0][1]*winding_currents[1] + \
                            L[0][2]*winding_currents[2] + L[0][3]*winding_currents[3]

    t1 += dt
