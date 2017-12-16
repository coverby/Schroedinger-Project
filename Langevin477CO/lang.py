import numpy as np
from .visualization import SimVis, start_server
import io
import asyncio

def read_energy(input_file):
    #Get the index, position, energy, and force from the position/energy file
    indx = []
    pos = []
    energy = []
    fx = []
    data = np.loadtxt(input_file)
    #print(data)
    #print(data[1][1])
    
    for row in range(len(data)):   
        indx.append(data[row][0])
        pos.append(data[row][1])
        energy.append(data[row][2])
        fx.append(data[row][3]) 
    return indx, pos, energy, fx


def read_coefficients(input_file):
    #Get the temperature, damping coefficient, tstep, and total integration time from file
    data = np.loadtxt(input_file)
    
    temp = data[0]
    damp = data[1]
    tstep = data[2]
    totaltime = data[3]
    return temp, damp, tstep, totaltime


def gdist(mean,temp,damp):
    #Generate Gaussian distributed eta for thermal disturbance
    sample = []
    var = np.multiply(np.multiply(2,temp),damp)
    for i in range(len(mean)):
        sample.append(np.random.normal(mean[i],var,1))
    return sample 

def write_output(idx, time, pos, vel, outname):
    #Write the output file containing index, time, position, and velocity of particles
    assert((len(idx) == len(time)) == (len(pos) == len(vel)))
    header_string = "#Index  Time    Position    Velocity"
    outfile = open(outname, "w")
    outfile.write(header_string)
    outfile.write("\n")
    for i in range(len(idx)):
        outfile.write("\t".join( (str(idx[i]), str(time[i]), str(pos[i]), str(vel[i])) ))
        outfile.write("\n")
    outfile.close()
    

def core_integrator(xi, vi, ui, fi, la, temp, mass, tstep, totaltime):
    #Core integrator using Velocity Verlet algorithm for Langevin equation
    totalsteps = int(np.floor(totaltime/tstep)) + 1
    assert(totalsteps > 2) #Check that we actually at least 2 time points
    xpos = np.zeros(totalsteps)
    accel = np.zeros(totalsteps)
    vel = np.zeros(totalsteps)
    potential = np.zeros(totalsteps)
    time = np.zeros(totalsteps)
    xpos[0] = xi
    accel[0] = fi/mass
    vel[0] = vi
    potential[0] = ui
    try: #An edge case here is if the position is exactly zero at the start
        kcoeff = potential[0]*2/(xpos[0]**2)
    except ZeroDivisionError:
        if np.isclose(0,potential[0]):
            kcoeff = 0 #Kind of a cludge, but garbage in garbage out
        else:
            kcoeff = potential[0]/(.001) #Also very much a cludge

    sample = gdist(time, temp, la)
    time[1] = tstep

    #Create timestep one for the Velocity Verlet Algorithm
    #x(t+dt) = x + v(t)*dt + 1/2*a(t)*dt^2
    #a(t+dt) from x(t+dt) using interaction potential
    #v(t+dt) = v(t) + (a(t) + a(t+dt))/2*dt

    xpos[1] = xpos[0] + vel[0]*tstep
    potential[1] = .5*kcoeff*xpos[1]**2
    accel[1] = (-la*vel[0] + sample[1]*tstep - 2*potential[1]/xpos[1])/mass
    vel[1] = vel[0] + (accel[0] + accel[1])*.5*tstep

    #Using the Velocity Verlet integration algorithm
    #Calculate new intermediate velocity, v(t+dt/2) = v + 1/2*a*dt
    #Calculate new position, x(t+dt) = x + v(t+dt/2)*dt
    #Calculate new acceleration from interaction potential
    #a(t+dt) = F(x(t+dt))
    #Calculate new end velocity
    #v(t+dt) = v(t+dt/2) + 1/2*a(t+dt)*dt

    for i in range(2,totalsteps):
        velint = vel[i-1] + .5*accel[i-1]*tstep
        xpos[i] = xpos[i-1] + velint*tstep

        potential[i] = .5*kcoeff*xpos[i]**2
        accel[i] = (-la*velint + sample[i]*tstep - kcoeff*xpos[i])/mass

        vel[i] = velint + .5*accel[i]*tstep
        time[i] = time[i-1] + tstep


    return xpos, accel, vel, potential, time

def main_handler(posfile, parafile, mass, vel, outfile, olength):
    #This function coordinates the other components after initialization
    idx, pos, energy, force = read_energy(posfile)
    temp, damp, tstep, totaltime = read_coefficients(parafile)

    idxout = []
    
    posout = []
    aclout = []
    velout = []
    potout = []
    timeout = []

    assert(olength<=int(np.floor(totaltime/tstep))), "Desired output length longer than expected time steps!"

    for i in range(len(idx)):
        output = core_integrator(pos[i],vel,energy[i],force[i],damp,temp,mass,tstep,totaltime)

        idxout.extend(np.multiply(np.ones(olength),idx[i]))
        posout.extend(output[0][-olength:])
        aclout.extend(output[1][-olength:])
        velout.extend(output[2][-olength:])
        potout.extend(output[3][-olength:])
        timeout.extend(output[4][-olength:])
        
    write_output(idxout,timeout,posout,velout,outfile)

def start(): #pragma: no cover
    sv = SimVis()
    print("This is a Langevin integrator utilizing the Verlet algorithm.")
    posfile = input("Enter the name of the file containing particle position information (format specified in readme): ")
    parafile = input("Enter the name of the file containing parameter information (format in readme): ")
    mass = float(input("Enter the particle mass: "))
    assert(mass > 0), "Mass must be a positive, non-zero value!"
    velocity = float(input("Enter particle velocity: "))
    outfile = input("Enter the name of the output file: ")
    olength = int(input("Enter how many timesteps for each particle you want in the output (enter 0 for all): "))
    main_handler(posfile, parafile, mass, velocity, outfile, olength)
