import numpy as np

def eclipseFractionAngle(rPlanet, rOrbit):
    return 2*np.arctan2(rPlanet, rOrbit)

def temperatureCoefficent(Voc, dVdT):
    """
    Solves for the temperature coefficient in terms of %/C
    """
    return (dVdT/Voc)*100

def truePowerFromRef(refPower, tempRef, tempCoef, tempActual):
    """
    Returns the true power of a solar cell given a reference power
    """
    return refPower*(1+tempCoef*(tempActual-tempRef))

def timeInEclipse(orbitalPeriod, eclipseFractionalAngle):
    return orbitalPeriod*(eclipseFractionalAngle/(2*np.pi))
    
def radiusOfOrbit(period):
    mu_e = 398600
    return ((mu_e*period**2)/(4*np.pi**2))**(1/3)

def powerOrbit(powerEclipse, timeEclipse, effEclipse, powerDay, timeDay, effDay):
    """
    Power required for one orbit
    """
    return (((powerEclipse*timeEclipse)/effEclipse)+(powerDay*timeDay)/effDay)/timeDay

# 10 year satellite mission
# in MEO
# period of 6 hours
# Panels normal to the sun
# Onboard Comp Power draw of 300W
# 50 extra watts during eclipse
# efficiency during day .8
# during night .9

def __main__():
    rOrbit = radiusOfOrbit(6*3600)
    print(rOrbit)
    pOrbit = powerOrbit(350, 41.66*60, .9, 300, 6*3600-41.66*60, .8)
    print(pOrbit)

if __name__ == "__main__":
    __main__()