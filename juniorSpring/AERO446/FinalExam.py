# ╔═════════════════════════════════════════════════════════════╗
# ║                                                             ║
# ║        .o.       oooooooooooo ooooooooo.     .oooooo.       ║
# ║       .888.      `888'     `8 `888   `Y88.  d8P'  `Y8b      ║
# ║      .8"888.      888          888   .d88' 888      888     ║
# ║     .8' `888.     888oooo8     888ooo88P'  888      888     ║
# ║    .88ooo8888.    888    "     888`88b.    888      888     ║
# ║   .8'     `888.   888       o  888  `88b.  `88b    d88'     ║
# ║  o88o     o8888o o888ooooood8 o888o  o888o  `Y8bood8P'      ║
# ║                                                             ║
# ║                ── CALIFORNIA POLYTECHNIC ──                 ║
# ║                                                             ║
# ╠═════════════════════════════════════════════════════════════╣
# ║   Author      :  Henry Flushman                             ║
# ║   Course      :  AERO446 - Spacecraft Electrical and        ║
# ║                                     Electric Systems        ║
# ║   Assignment  :  Final Exam                                 ║
# ║   Date        :  June 8, 2026                               ║
# ╚═════════════════════════════════════════════════════════════╝


# === imports ===
import numpy as np

# === constants ===
Re = 6378.18e3
c = 3.0e8


# === Helper Functions ===

def linear_to_db(x):
    return 10*np.log10(x)

def db_to_linear(x):
    return 10**(x/10)

def frequency_to_wavelength(f):
    return c / f

def characteristic_dimensions(h, Re=Re):
    R = h + Re
    FOR = 2*np.arcsin(Re/R)
    d_edge = Re/np.tan(0.5*FOR)
    y = d_edge*np.sin(0.5*FOR)
    z = np.sqrt(Re**2 - y**2)
    x = R - z
    return x, y, z, d_edge, np.rad2deg(FOR)

def transmitterDiameter(wavelength, FOR_deg):
    return (70*wavelength)/FOR_deg

def areaFromDiameter(diameter):
    return (np.pi*diameter**2)/4

def transmitterGain(transmitterArea, wavelength, eta=0.55):
    return (4*np.pi*eta*transmitterArea)/wavelength**2

def freeSpacePathLoss(pathLength, wavelength):
    return -20*np.log10((4*np.pi*pathLength)/wavelength)

def freeSpacePathLoss2(pathlength, wavelength):
    return 20*np.log10(wavelength/(4*np.pi*pathlength))

def atmosphericLoss(loss_per_km, phi_min):
    angle = np.deg2rad(90-phi_min)
    return 100/np.cos(angle) * loss_per_km

def systemTemperatureGain(noiseFactor, tempAntenna=290):
    return (noiseFactor-1)*290 + tempAntenna

def recieverArea(antennaGain, wavelength, eta=0.55):
    return (antennaGain*wavelength**2)/(4*np.pi*eta)
    


# === Problem 1 ===


# given

R_geo = 42164.18e3  # m
f = 40e9            # Hz
B = 10e6            # Hz
P_tx = 60.0         # W
OBO = 1.0           # dB
F_rx = 2.0          # dB
eta = 0.55          
L_o2_h2o = 0.065    # dB/K
C_over_N = 3        # dB

wavelength = frequency_to_wavelength(f)
P_tx = linear_to_db(P_tx)
F_rx = db_to_linear(F_rx)

# assumed
L_polar = 0.3
L_radome = 1.0


# Analysis

print("="*50)
print("- Problem 1")
print("="*50)

x, y, z, d_edge, FOR_deg = characteristic_dimensions(h=R_geo-Re)     

D_tx = transmitterDiameter(wavelength, FOR_deg)
A_tx = areaFromDiameter(D_tx)
G_tx = transmitterGain(A_tx, wavelength)

G_tx = linear_to_db(G_tx)

EIRP = P_tx + G_tx - OBO

print("1.1: EIRP\n")
print(f"     {EIRP:.3f} dB\n")

L_fs = freeSpacePathLoss(d_edge, wavelength)

L_atm = atmosphericLoss(0.065, 5)

L_tot = L_fs - L_atm - L_polar - L_radome

print("1.2: Total Path Loss\n")
print(f"     {L_tot:.3f} dB\n")

T_sys = systemTemperatureGain(F_rx)

print("1.3: Ground System Noise Temperature\n")
print(f"     {T_sys:.3f} dB-K\n")

G_over_T = C_over_N - EIRP - L_tot - 228.6 + linear_to_db(B)

G_rx = G_over_T + linear_to_db(T_sys)

print("1.4: Total Ground System Gain for C/N = 3dB\n")
print(f"     {G_rx:.3f} dB\n")

G_antenna = db_to_linear(G_rx - 20)

A_rx = recieverArea(G_antenna, wavelength)

D_rx = 2*np.sqrt(A_rx/np.pi)

print("1.5: Diameter of Receiving Antenna\n")
print(f"     {D_rx:.3f} m\n")

print("="*50)

# === Problem 2 ===

# 2.1
# In general, as frequency increases atmospheric attenuation decreases

# False, atmospheric attentuation is dependent on frequency. As frequency increases, the wall penetrative
# abilities of electromagnetic waves decrease as compared to lower
# frequencies which can travel further unobstructed.


# 2.2
# For parabolic antennas, as frequency increases free space loss increases

# True, as frequency increases you will incur more free space path losses. 
# This can be proved by the nature of the equation, having a negative sign 
# infront of log with wavelength in the denominator implies that an increase
# in wavelength decreases the magnitude of the equation


# 2.3
# The K-band is used for defense applications

# False, while the K-band is still sometimes used by the defense industry
# it is the X-band and Ka-band that are most commonly used


# 2.4
# The noise temperature of satellite’s transmitting antenna is affected by the field of view

# False, the receiving antenna noise temperature is affected by what the antenna
# can see in its field of view. A transmitting antenna doesn't really have a noise
# temperature in a link-budget sense.


# 2.5
# The gain of an antenna communicating in the S-band is higher than the gain
# of the same antenna communicating in the X-band

# False, for the same physical antenna, gain will increase with frequency.
# Because the X-band has a higher frequency than S-band, the antenna gain
# is higher in the X-band


# === Problem 3 ===

print("- Problem 3")
print("="*50)

def hamming_distance(a, b):
    return sum(bit1 != bit2 for bit1, bit2 in zip(a, b))

codebook = {
    "00": "00000",
    "01": "00101",
    "10": "10111",
    "11": "01111"
}

received_message = [
    "00101",
    "10101",
    "10111",
    "11111",
    "01111",
    "00001"
]

decoded_options = []

print("Received Word | Distances to Codebook        | Decoded Data")
print("-------------------------------------------------------------")

for received in received_message:
    distances = {}

    for data, codeword in codebook.items():
        distances[data] = hamming_distance(received, codeword)

    min_distance = min(distances.values())

    closest_data = [
        data for data, distance in distances.items()
        if distance == min_distance
    ]

    decoded_options.append(closest_data)

    distance_list = [distances[data] for data in codebook.keys()]

    if len(closest_data) == 1:
        decoded_text = closest_data[0]
    else:
        decoded_text = " or ".join(closest_data)

    print(f"{received:13} | {distance_list}                 | {decoded_text}")
    
    
possible_messages = [""]

for options in decoded_options:
    new_messages = []

    for message in possible_messages:
        for option in options:
            new_messages.append(message + option)

    possible_messages = new_messages


print("\nDecoded pieces:")
for received, options in zip(received_message, decoded_options):
    print(f"{received} -> {' or '.join(options)}")

print("\nPossible original transmitted messages:")
for message in possible_messages:
    print(message)