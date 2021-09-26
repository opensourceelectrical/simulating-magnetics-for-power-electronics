import math

dt = 1.0e-6

r = 0.01    # Winding resistance

# Coil/core construction
no_of_turns = 200
cs_area = 9.0e-4
length_iron = 15.0e-2
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
    R_iron = length_iron / (mu * cs_area)
    # phi = NI / R
    # I = phi * R / N
    ind_current = flux * R_iron / no_of_turns

    vsrc = vmeas - ind_current*series_res

    indmodel_flux = flux
    indmodel_emf = vmeas - ind_current*r

    t1 += dt
