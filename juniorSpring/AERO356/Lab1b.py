import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Load CSV
file_path = 'juniorSpring/AERO356/labData/Lab1b_radiation.csv'
df = pd.read_csv(file_path)

runs = []
num_runs = 5
time_offset = 0.0
run_boundaries = []

# -------------------------------
# BUILD df_all (THIS FIXES YOUR ERROR)
# -------------------------------
for r in range(1, num_runs + 1):

    time_col = f'Time (s) Run #{r}'
    temp1_col = f'Temperature 1, Ch A:1 (°C) Run #{r}'
    temp2_col = f'Temperature 2, Ch A:1 (°C) Run #{r}'
    temp3_col = f'Temperature 3, Ch A:1 (°C) Run #{r}'
    temp4_col = f'Temperature 2, Ch B:1 (°C) Run #{r}'
    temp5_col = f'Temperature 3, Ch B:1 (°C) Run #{r}'

    required_cols = [time_col, temp1_col, temp2_col, temp3_col, temp4_col, temp5_col]

    if not all(col in df.columns for col in required_cols):
        continue

    run_df = df[required_cols].copy()
    run_df.columns = [
        'Time (s)',
        'Temp1_ChA',
        'Temp2_ChA',
        'Temp3_ChA',
        'Temp2_ChB',
        'Temp3_ChB'
    ]

    run_df = run_df.dropna().reset_index(drop=True)

    if run_df.empty:
        continue

    run_df['Time (s)'] += time_offset
    runs.append(run_df)

    dt = run_df['Time (s)'].iloc[1] - run_df['Time (s)'].iloc[0]
    run_boundaries.append(run_df['Time (s)'].iloc[-1])
    time_offset = run_df['Time (s)'].iloc[-1] + dt

df_all = pd.concat(runs, ignore_index=True)

# -------------------------------
# PLOT
# -------------------------------
plt.figure(figsize=(12, 7))

plt.plot(df_all['Time (s)'], df_all['Temp1_ChA'], label='Large Black Aluminum')
plt.plot(df_all['Time (s)'], df_all['Temp2_ChA'], label='Black Aluminum')
plt.plot(df_all['Time (s)'], df_all['Temp3_ChA'], label='Black Brass')
plt.plot(df_all['Time (s)'], df_all['Temp2_ChB'], label='Polished Aluminum')
plt.plot(df_all['Time (s)'], df_all['Temp3_ChB'], label='White Aluminum')

for boundary in run_boundaries[:-1]:
    plt.axvline(boundary, linestyle='--')

plt.xlabel('Time (s)')
plt.ylabel('Temperature (°C)')
plt.title('Cylinder Temperature vs Time')
plt.legend()
plt.grid(True)
plt.show()

# -------------------------------
# CALCULATIONS
# -------------------------------

# Specific heat [J/kg-K]
specific_heat = {
    'Large Black Aluminum': 900,
    'Black Aluminum': 900,
    'Black Brass': 380,
    'Polished Aluminum': 900,
    'White Aluminum': 900
}

# Convert units properly:
# grams → kg, inches → meters
inch_to_m = 0.0254

properties = {
    'Large Black Aluminum': {'mass': 52.690/1000, 'length': 1.5*inch_to_m, 'diameter': 1.0*inch_to_m},
    'Black Aluminum': {'mass': 29.510/1000, 'length': 1.5*inch_to_m, 'diameter': 0.75*inch_to_m},
    'Black Brass': {'mass': 91.498/1000, 'length': 1.5*inch_to_m, 'diameter': 0.75*inch_to_m},
    'Polished Aluminum': {'mass': 27.741/1000, 'length': 1.5*inch_to_m, 'diameter': 0.75*inch_to_m},
    'White Aluminum': {'mass': 29.569/1000, 'length': 1.5*inch_to_m, 'diameter': 0.75*inch_to_m}
}

column_map = {
    'Large Black Aluminum': 'Temp1_ChA',
    'Black Aluminum': 'Temp2_ChA',
    'Black Brass': 'Temp3_ChA',
    'Polished Aluminum': 'Temp2_ChB',
    'White Aluminum': 'Temp3_ChB'
}

# Choose EARLY heating region (IMPORTANT)
fit_start = 0
fit_end = 300   # you can tweak this

results = []

for cylinder, temp_col in column_map.items():

    fit_data = df_all[(df_all['Time (s)'] >= fit_start) &
                      (df_all['Time (s)'] <= fit_end)]

    time = fit_data['Time (s)'].values
    temp = fit_data[temp_col].values

    slope, _ = np.polyfit(time, temp, 1)

    m = properties[cylinder]['mass']
    L = properties[cylinder]['length']
    D = properties[cylinder]['diameter']
    A = L * D
    c = specific_heat[cylinder]

    I = (m * c / A) * slope

    results.append([cylinder, slope, m, A, I])

results_df = pd.DataFrame(results, columns=[
    'Cylinder',
    'Slope dT/dt [C/s]',
    'Mass [kg]',
    'Projected Area [m^2]',
    'Intensity [W/m^2]'
])

print(results_df)