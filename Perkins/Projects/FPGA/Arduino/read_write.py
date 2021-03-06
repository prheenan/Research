# force floating point division. Can still use integer with //
from __future__ import division
# other good compatibility recquirements for python3
from __future__ import absolute_import
from __future__ import print_function
from __future__ import unicode_literals
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys


import serial
import time

def run():
    """
    shamelessly copied from:
    
    petrimaki.com/2013/04/28/reading-arduino-serial-ports-in-windows-7/
    """
    port_str = "/dev/cu.usbmodem1411"
    ser = serial.Serial(port_str, 9600, timeout=0)
    # sleep, to allow for initialization (see: 
    # https://playground.arduino.cc/Interfacing/Python
    # especially:
    # "work around this issue by simply placing a 'time.sleep(2)' call 
    # between the serial connection and the write call."
    time.sleep(3)
    # forever read/write 
    while 1:
        str_input = raw_input("Enter a character: ")
        float_input = float(str_input)
        sstr='{0}\n'.format(float_input)
        print("sending {:s}".format(sstr))
        ser.write(sstr.encode('ascii'))
        try:
            time.sleep(2)
            line = ser.readline()
            # XXX why is this necessary? arduino writes an extra line?
            _ = ser.readline()
        except ser.SerialTimeoutException:
            print('Data could not be read')
        print("Arduino sent back:\n\t[{:s}]".format(line.strip()))

if __name__ == "__main__":
    run()
