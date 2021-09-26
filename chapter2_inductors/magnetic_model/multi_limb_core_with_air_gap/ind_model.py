import math

dt = 1.0e-6

r = 0.01        # winding resistance

# Coil/core construction
no_of_turns = 200
cs_area = 9.0e-4
length_a = 5.0e-2
length_b = 4.5e-2
length_c = 4.0e-2
length_w = 3.0e-2
length_lg_2 = 0.1e-3
mu_0 = 4*math.pi*1.0e-7
mu_r = 1000.0
mu = mu_0 * mu_r

if t_clock >= t1:

    # dpsi/dt = v - ir
    # Runge Kutta Fourth Order method
    k1 = (vmeas - ind_current*r)
    k2 = (vmeas - (ind_current + dt*k1/2.0)*r)
    k3 = (vmeas - (ind_current + dt*k2/2.0)*r)
    k4 = (vmeas - (ind_current + dt*k3)*r)
    k = (k1 + k2*2 + k3*2 + k4)*dt/6.0
    flux_linkage += k

    # phi = psi / N
    flux = flux_linkage / no_of_turns
    # R = l / (mu A)
    R1 = (2*length_a + length_b + 2.0*length_w \
            - 2*length_lg_2) / (mu * cs_area)
    Rg1 = 2*length_lg_2 / (mu_0 * cs_area)
    R2 = (length_b + length_w) / (mu * cs_area)
    R3 = (2*length_c + length_b + 2.0*length_w \
            - 2*length_lg_2) / (mu * cs_area)
    Rg2 = length_lg_2 / (mu_0 * cs_area)
    R_eq = R1 + Rg1 + (R2 * (R3 + Rg2) / (R2 + R3 + Rg2))
    # phi = NI / R
    # I = phi * R / N
    ind_current = flux * R_eq / no_of_turns

    vsrc = vmeas - ind_current*series_res

    indmodel_emf = vmeas - ind_current*r
    indmodel_flux = flux
    indmodel_flux2 = flux * (R2 * (R3 + Rg2) / (R2 + R3 + Rg2)) / R2
    indmodel_flux3 = flux * (R2 * (R3 + Rg2) / (R2 + R3 + Rg2)) / R3

    t1 += dt
