import serial

def open_port(port_name):
    print(port_name)
    ser = serial.Serial(port_name,9600)
    print("Opened port "+port_name)
    port_ok = True

    if (not(port_ok)):
        raise Exception("Failed to open port")
    return ser