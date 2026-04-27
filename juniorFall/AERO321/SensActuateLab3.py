import os
import serial as ser
import time as systime
import math
import ArduinoFunc as ard
import numpy as np
import scipy.stats
from scipy.stats import chi2
import matplotlib.pyplot as plt

""" Sensors and Actuators Lab 2 """

ser = ard.open_port('COM11')

data = []
while len(data) < 200:  # read ~10s at 10Hz
    line = ser.readline().decode().strip().split(',')
    data.append([float(x) for x in line])
ser.close()

data = np.array(data)
dist = data[:,1]
motorspeed = data[:,2]

print(f"Distance: {dist} \nSpeed: {motorspeed}")

time = np.linspace(1,20,200)

plt.figure()
plt.plot(time, dist)
plt.xlabel('Distance (m)')
plt.ylabel('Time (s)')
plt.title('Distance vs Time')

plt.figure()
plt.plot(time, motorspeed)
plt.xlabel('Speed (rad/s)')
plt.ylabel('Time (s)')
plt.title('Speed vs Time')
plt.show()