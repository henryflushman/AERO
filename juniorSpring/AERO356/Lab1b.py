import pandas as pd
import matplotlib.pyplot as plt

# Load the CSV file
file_path = 'AERO356/labData/Lab1b_radiation.csv'
df = pd.read_csv(file_path)

runs = []
num_runs = 5
time_offset = 0.0
run_boundaries = []

for r in range(1, num_runs + 1):
    # Column names for this run
    time_col = f'Time (s) Run #{r}'
    temp1_col = f'Temperature 1, Ch A:1 (°C) Run #{r}'
    temp2_col = f'Temperature 2, Ch A:1 (°C) Run #{r}'
    temp3_col = f'Temperature 3, Ch A:1 (°C) Run #{r}'
    temp4_col = f'Temperature 2, Ch B:1 (°C) Run #{r}'
    temp5_col = f'Temperature 3, Ch B:1 (°C) Run #{r}'

    required_cols = [time_col, temp1_col, temp2_col, temp3_col, temp4_col, temp5_col]

    # Skip missing runs safely
    if not all(col in df.columns for col in required_cols):
        print(f'Run #{r} not found, skipping.')
        continue

    # Extract and rename to common names
    run_df = df[required_cols].copy()
    run_df.columns = [
        'Time (s)',
        'Temp1_ChA',
        'Temp2_ChA',
        'Temp3_ChA',
        'Temp2_ChB',
        'Temp3_ChB'
    ]

    # Drop rows with missing data
    run_df = run_df.dropna().reset_index(drop=True)

    # Skip empty runs
    if run_df.empty:
        print(f'Run #{r} is empty after cleaning, skipping.')
        continue

    # Append in time
    run_df['Time (s)'] = run_df['Time (s)'] + time_offset
    run_df['Run'] = r
    runs.append(run_df)

    # Estimate time step for clean continuation
    if len(run_df) > 1:
        dt = run_df['Time (s)'].iloc[1] - run_df['Time (s)'].iloc[0]
    else:
        dt = 0.0

    run_boundaries.append(run_df['Time (s)'].iloc[-1])
    time_offset = run_df['Time (s)'].iloc[-1] + dt

# Combine all runs
df_all = pd.concat(runs, ignore_index=True)

# Plot appended data
plt.figure(figsize=(12, 7))

plt.plot(df_all['Time (s)'], df_all['Temp1_ChA'], label='Temperature 1, Ch A:1', linewidth=2)
plt.plot(df_all['Time (s)'], df_all['Temp2_ChA'], label='Temperature 2, Ch A:1', linewidth=2)
plt.plot(df_all['Time (s)'], df_all['Temp3_ChA'], label='Temperature 3, Ch A:1', linewidth=2)
plt.plot(df_all['Time (s)'], df_all['Temp2_ChB'], label='Temperature 2, Ch B:1', linewidth=2)
plt.plot(df_all['Time (s)'], df_all['Temp3_ChB'], label='Temperature 3, Ch B:1', linewidth=2)

# Add vertical lines between runs
for boundary in run_boundaries[:-1]:
    plt.axvline(boundary, linestyle='--', linewidth=1, color='black')

plt.xlabel('Time (s)')
plt.ylabel('Temperature (°C)')
plt.title('Temperature vs Time (Appended Runs)')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()