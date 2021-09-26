dt = 1.0e-6

L = 0.3  # Inductance of coil
r = 0.01 # Winding resistance

if t_clock >= t1:

    # Backward Euler method
    # ind_current += (1/L) * (vmeas - ind_current*r)*dt
    # Runge Kutta Fourth Order method
    k1 = (1/L) * (vmeas - ind_current*r)
    k2 = (1/L) * (vmeas - (ind_current + dt*k1/2.0)*r)
    k3 = (1/L) * (vmeas - (ind_current + dt*k2/2.0)*r)
    k4 = (1/L) * (vmeas - (ind_current + dt*k3)*r)
    k = (k1 + k2*2 + k3*2 + k4)*dt/6.0
    ind_current += k

    vsrc = vmeas - ind_current*series_res

    inductor_emf = vmeas - ind_current * r
    inductor_fluxlinkage = ind_current * L

    t1 += dt
