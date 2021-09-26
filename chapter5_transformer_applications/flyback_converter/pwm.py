import math

dt = 1.0e-9
sw_freq = 10000.0  # 10 kHz
duty_ratio = 0.25

if t_clock >= t1:

    # starts at 0 and reaches 1 every sw cycle
    # changes from 0 to 1 in (1/sw_freq)
    # Slope = 1 / (1/sw_freq) = sw_freq
    carr_wave += (sw_freq)*dt

    if carr_wave > 1.0:
        carr_wave = 0.0

    if duty_ratio > carr_wave:
        sw_control = 1.0
    else:
        sw_control = 0.0

    t1 += dt

pwm_carr_wave = carr_wave
pwm_duty_ratio = duty_ratio
pwm_sw_control = sw_control
