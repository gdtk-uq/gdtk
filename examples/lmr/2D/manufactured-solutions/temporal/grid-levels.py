gridLevels = [
    {'dt': 1.0e-5},
    {'dt': 0.5e-5},
    {'dt': 0.25e-5},
    {'dt': 0.125e-5},
    {'dt': 0.0625e-5},
]

for gl in gridLevels:
    gl['dx'] = gl['dt']

