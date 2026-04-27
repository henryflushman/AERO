#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Simple example of reading serial data from an Arduino
"""

import serial
import statistics
import numpy as np
import matplotlib.pyplot as plt
 
# try to find a serial port, this may need to be changed on Windows
# make sure to close Serial Monitor/Plotter in Arduino
port_ok = False
for i in range(5,6):
    try:
        port_name = "COM%d"%i
        #alternatives /dev/ttyUSB%d, COM%d, /dev/ttyS%d
        
        print(port_name)
        ser = serial.Serial(port_name,9600)
        print("Opened port "+port_name)
        port_ok = True
        break
    except serial.SerialException:
        pass
 
if (not(port_ok)):
    raise Exception("Failed to open port")
    
vals = [0]
res = [0]
temp = [0]
for i in range(150):
    line = str(ser.readline()); # read line from serial port
    line = line[2:-5]  #eliminate trailing b' and \r\n
    pieces = line.split(','); # split by comma
    if i == 0:
        vals = [float(pieces[0])]
        res = [float(pieces[1])]
        temp = [float(pieces[2])]
    else:
        vals.append(float(pieces[0]))
        res.append(float(pieces[1]))
        temp.append(float(pieces[2]))


    print(pieces)
    
print(temp)
prcn = (max(vals) - min(vals))/1023.0 * 100.0
print("Range: %g, which is %g percent"%(max(vals)-min(vals),prcn))
avg_temp = statistics.mean(temp)
temp_acc = (avg_temp - 22.3)/22.3 * 100.0
print(f"Average temperature: {avg_temp:.5f} C, which is {abs(temp_acc):.5f} percent from 22.3 C")

span = max(temp) - min(temp)
firstorder_temp = 0.632*span + min(temp)

print(f"Range of temperature: {[min(temp), max(temp)]}")

ser.close()

time = np.linspace(1,30,150)

import numpy as np

t = np.asarray(time, dtype=float)
y = np.asarray(temp, dtype=float)
Tstar = float(firstorder_temp)

if Tstar < y.min() or Tstar > y.max():
    print("firstorder_temp is outside the data range; no crossing.")
else:
    below = y < Tstar
    crossing_idxs = np.where(below[:-1] & ~below[1:])[0]

    if crossing_idxs.size:
        k = int(crossing_idxs[0]) 
        t0, t1 = t[k], t[k+1]
        y0, y1 = y[k], y[k+1]
        t_cross = t0 + (Tstar - y0) * (t1 - t0) / (y1 - y0)
        print(f"First-order level at t ≈ {t_cross:.2f} s (between samples {k} and {k+1})")
    else:
        print("No crossing found (data may never cross or sits flat at threshold).")


plt.plot(time, temp)
plt.axhline(Tstar, linestyle='--')
plt.axvline(t_cross, linestyle=':')
plt.xlabel("Time [seconds]")
plt.ylabel("Temperature [C]")
plt.grid(True)
plt.show()


input("Type anything to continue:")

import os

os.system('cls')



""" Plot and analyze the range of a Thermistor """

ADC = np.linspace(5,1018,1013)
R = 10000/((1023/ADC) - 1)
T = (1/298.15 + (1/3380)*np.log(R/10000))**-1 - 273.15

plt.scatter(ADC, T, s=1)
plt.grid(True)
plt.xlabel("ADC Value [5:1018]")
plt.ylabel("Temperature [C]")
plt.show()