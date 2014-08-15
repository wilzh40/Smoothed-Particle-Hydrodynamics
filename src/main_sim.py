import math
import os
import numpy as np
import sys
from pylab import *

def parallel_waves(n=24, #26 for our first test?
                   time=0, 
                   phi=math.pi,
                   amplitude=1,
                   velocity=0.0001):
    """
    Array of two travelling waves, second one starts
    half way through the array
    """

    if n % 2 != 0:
        raise NotImplementedError("Currently only supports even number of muscles!")

    j = n/2

    row_positions = np.linspace(0,1.5*2*math.pi,j)

    wave_1 = (map(math.sin,(row_positions - velocity*time)))
    wave_2 = (map(math.sin,(row_positions + (math.pi) - velocity*time)))

    normalize_sine = lambda x : (x + 1)/2
    wave_1 = map(normalize_sine, wave_1)
    wave_2 = map(normalize_sine, wave_2)

    double_wave_1 = []
    double_wave_2 = []

    for i in wave_1:
        double_wave_1.append(i)
        double_wave_1.append(i)

    for i in wave_2:
        double_wave_2.append(i)
        double_wave_2.append(i)
        
    return (double_wave_1,double_wave_2)
    
def input_neural_data():
	volts = {}
	traces = open('c302_A.dat','r')

	print ("Loaded .dat file")
	# Very inefficient...
	for line in traces:
		if not line.strip().startswith('#'):
			points = line.split()
			for i in range(len(points)):
				if not volts.has_key(i):
					volts[i] = []
				volts[i].append(float(points[i])+0.0020*i)

	
	return volts
        
class muscle_simulation():

    def __init__(self,increment=10.0):
		os.getcwd()
		self.increment = increment
		self.t = 0
		self.signal_array = input_neural_data()
   
        

    def run(self,do_plot = True):
        self.contraction_array =  parallel_waves(time = self.t)
        self.t += self.increment
        return list(np.concatenate([self.contraction_array[0],
                                    self.contraction_array[1],
                                    self.contraction_array[1],
                                    self.contraction_array[0]]))
        #return(self.contraction_array)
        
