dt = 1.0e-8         # Control time step

f_sw = 20000.0      # 20 kHz
T_sw = 1 / f_sw     # Switching time period
carr_slope = 1 / T_sw

if t_clock >= t1:
    # Carrier (sawtooth) waveform
    carr_wave += carr_slope*dt
    if carr_wave > 1.0:
        carr_wave = 0.0

    duty_ratio = 0.3
    if duty_ratio > carr_wave:
        s1gate = 1.0
    else:
        s1gate = 0.0

    pwm_carr = carr_wave
    pwm_dutyratio = duty_ratio
    pwm_gate = s1gate

    t1 += dt
