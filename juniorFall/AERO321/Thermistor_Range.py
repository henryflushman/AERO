import os
os.system('cls')

import numpy as np
import matplotlib.pyplot as plt

""" Plot and analyze the range of a Thermistor """

ADC = np.linspace(5,1018,1013)
R = 10000/((1023/ADC) - 1)
T = (1/298.15 + (1/3380)*np.log(R/10000))**-1 - 273.15

plt.scatter(ADC, T, s=1)
plt.grid(True)
plt.xlabel("ADC Value [5:1018]")
plt.ylabel("Temperature [C]")
plt.show()