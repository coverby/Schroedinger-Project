import io
import numpy as np
import Langevin477CO as lang
import os

def test_read_energy_read():
    '''Tests if the function reads and returns given an input file'''
    test_string = '''
    #test input energy
    #i  x   energy  fx 
    0   0   0   -2
    1   0   1   -1
    2   0   2   4
    3   0   4   6
    4   3   0   1
    '''

    test_file = io.StringIO(test_string)
    indx, pos, energy, fx = lang.read_energy(test_file)
    #print(pos)
    #print(energy)
    #print((np.isclose(pos, [0,1,2,3,4])).any())
    assert((np.isclose(indx, [0,1,2,3,4])).any())
    assert((np.isclose(pos, [0,0,0,0,3])).any())
    assert((np.isclose(energy,[0,1,2,4,0])).any())
    assert((np.isclose(fx, [-2,-1,4,6,1])).any())


def test_read_coefficients_read():
    test_string = '''
    #Test coefficient string
    #temp   damping    tstep    totaltime 
    273 .123    0.001   10
    '''
    test_file = io.StringIO(test_string)
    temp, damp, tstep, totaltime = lang.read_coefficients(test_file)

    assert(np.isclose(temp,273))
    assert(np.isclose(damp, .123))
    assert(np.isclose(tstep,0.001))
    assert(np.isclose(totaltime, 10))

def test_write_output():
    #Test the ability to write correctly formatted output text files.
    test_data = [[1,10.0,-.432,.234],[2,10.0,123,456]]
    test_idx = [1, 2]
    test_time = [10.0, 10.0]
    test_pos = [-.432, 123]
    test_vel = [.234, 456]
    fname = "writetest.txt"
    lang.write_output(test_idx,test_time,test_pos,test_vel,fname)

    file_data = np.loadtxt(fname)
    assert(np.isclose(file_data, test_data).any())
    os.remove(fname)

def test_main_handler():
    #Test the main handler function's ability to coordinate output file format
    posstring = '''
    #test input energy
    #i  x   energy  fx 
    0   0.1   0.1   -2
    1   0.1   1   -1
    2   0.2   2   4
    3   -2   4   6
    4   3   0   1
    '''
    parastring = '''
    #Test coefficient string
    #temp   damping    tstep    totaltime 
    273 .123    0.01   10
    '''
    mass = 1
    vel = 1
    olength = 3
    posfile = io.StringIO(posstring)
    parafile = io.StringIO(parastring)
    outfile = "handlertestout.txt"
    lang.main_handler(posfile, parafile, mass, vel, outfile, olength)
    
    file_data = np.loadtxt(outfile)
    assert(len(file_data) == 15)
    assert(file_data[14][0] == 4)
    os.remove(outfile)

def test_temp_distribution():
    #Test to make sure that we have a distribution with the right mean
    nsam = 1000
    temp = 300
    damp = .123
    mean = np.zeros(nsam)
    samples = lang.gdist(mean, temp, damp)
    print(len(samples))
    if np.absolute(np.average(samples)) >= 5:
        newavg = 0
        for i in range(3):
            samples1 = lang.gdist(mean, temp, damp)
            newavg += np.average(samples1)
            samples1 = []
        assert(newavg <= 15)

def test_core_integrator_setup():
    #Test that the integrator can be loaded and outputs in the expected format
    xi = 1
    vi = 1
    ui = 1
    fi = 1
    la = .123
    temp = 300
    mass = 1
    tstep = .01
    totaltime = 10
    output = lang.core_integrator(xi,vi,ui,fi,la,temp,mass,tstep,totaltime)
    print("Number of output variables: ")
    print(len(output))
    assert(len(output) == 5)
    print("Length of each output: ")
    print(len(output[0]))
    assert(len(output[0]) == int(np.floor(totaltime/tstep)+1 ))

def test_integrator_conserve_potential():
    #Solvent interactions + potential = conserves expected potential energy
    #With no energy input from random fluctuations, we expect total energy to decay to 0
    xi = 1
    vi = 1
    ui = 1
    fi = 1
    la = 0.7
    temp = 0 #Basically disables the random inputs (var => 0)
    mass = 1
    tstep = .001
    totaltime = 200
    output = lang.core_integrator(xi,vi,ui,fi,la,temp,mass,tstep,totaltime)
    #output is xpos, accel, velocity, potential energy, and time
    print("Potential energy format looks like:")
    print(output[3])
    print("Average potential energy over time looks like:")
    print(np.average(output[3]))
    assert(np.isclose(np.average(output[3]),0,atol=1.e-2))

def test_integrator_conserve_KE():
    #Considering only potential = conserved KE
    #Without any solvent interactions or thermal input, the integrator must conserve total energy
    xi = 1
    vi = 1
    ui = 1
    fi = 1
    la = 0 #Disables both solvent interactions (no drag) and random inputs (var => 0)
    temp = 300
    mass = 1
    tstep = .001
    totaltime = 200
    output = lang.core_integrator(xi,vi,ui,fi,la,temp,mass,tstep,totaltime)
    #output is xpos, accel, velocity, potential energy, and time
    # avgKE = .5*mass*vel^2
    sum = 0
    for i in range(len(output[2])):
        sum += .5*mass*output[2][i]**2
    avgKE = sum/len(output[2])
    print("Initial KE:")
    print(ui+.5*mass*vi**2)
    print("Average kinetic energy over time looks like:")
    print(np.average(output[3]) + avgKE)

    assert(np.isclose(np.average(output[3] + avgKE),(ui+.5*mass*vi**2),rtol=1.e-2))