import pandas as pd
import matplotlib.pyplot as plt

# Load the CSV file
file_path = 'AERO356/labData/Lab1a_bottles.csv'
df = pd.read_csv(file_path)

# Keep only the columns we need
cols = [
    'Time (s) Run #2',
    'Black Bottle',
    'Brushed Aluminum Bottle',
    'White Bottle'
]
df = df[cols].copy()

# Remove rows with missing temperature values.
# In this dataset, that effectively skips every other row,
# while also being safer in case the missing-value pattern changes.
df = df.dropna(subset=['Black Bottle', 'Brushed Aluminum Bottle', 'White Bottle'])

# Optional: reset row numbers after cleaning
df = df.reset_index(drop=True)

# Extract data for plotting
time = df['Time (s) Run #2']
black = df['Black Bottle']
aluminum = df['Brushed Aluminum Bottle']
white = df['White Bottle']

# Create the plot
plt.figure(figsize=(10, 6))
plt.plot(time, black, label='Black Bottle', linewidth=2)
plt.plot(time, aluminum, label='Brushed Aluminum Bottle', linewidth=2)
plt.plot(time, white, label='White Bottle', linewidth=2)

# Labels, title, legend, and grid
plt.xlabel('Time (s)')
plt.ylabel('Temperature (°C)')
plt.title('Water Bottle Temperature vs Time')
plt.legend()
plt.grid(True)
plt.tight_layout()

# Show the plot
plt.show()