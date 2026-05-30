import pandas as pd
import matplotlib.pyplot as plt

file_path = 'juniorSpring/AERO356/labData/Lab1a_bottles.csv'
df = pd.read_csv(file_path)

time = df['Time (s) Run #2']

plt.figure(figsize=(10, 6))

# Black bottle
plt.plot(time[::40], df['Black Bottle'][::40],
         linestyle='None', marker='o',
         markersize=4, color='black',
         label="Black Bottle")

# Aluminum bottle
plt.plot(time[::40], df['Brushed Aluminum Bottle'][::40],
         linestyle='None', marker='s',
         markersize=4, color='black',
         label="Brushed Aluminum Bottle")

# White bottle
plt.plot(time[::40], df['White Bottle'][::40],
         linestyle='None', marker='^',
         markersize=4, color='black',
         label="White Bottle")

plt.xlabel('Time (s)')
plt.ylabel('Temperature (°C)')
plt.title('Water Bottle Temperature vs Time')

plt.legend()
plt.grid(True, linestyle=':', linewidth=0.5)
plt.tight_layout()
plt.show()