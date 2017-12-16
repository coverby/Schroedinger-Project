import numpy as np
from .visualization import SimVis, start_server
import io
import asyncio

def read_energy(input_file):
    #datafile = io.StringIO(input_file)
    #print(input_file.getvalue())
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
    data = np.loadtxt(input_file)
    
    temp = data[0]
    damp = data[1]
    tstep = data[2]
    totaltime = data[3]
    return temp, damp, tstep, totaltime


def gdist(mean,temp,damp):
    #Generate Gaussian distribuuted nu
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
    totalsteps = int(np.floor(totaltime/tstep))
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
    kcoeff = potential[0]*2/(xpos[0]**2)
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

def main_handler(posfile, parafile, mass, outfile, olength):
    #This function coordinates the other components after initialization

def start():  #pragma: no cover
    print("This is a Langevin integrator utilizing the Verlet algorithm.")
    posfile = input("Enter the name of the file containing particle position information: ")
    parafile = input("Enter the name of the file containing parameter information: ")
    mass = input("Enter the particle mass: ")
    assert(mass > 0), "Mass must be a positive, non-zero value!"
    outfile = input("Enter the name of the output file: ")
    olength = input("Enter how many timesteps for each particle you want in the output (enter 0 for all): ")
    assert(int(olength)), "Enter an integer for output length"
    main_handler(posfile, parafile, mass, outfile, olength)

#def start(): #pragma: no cover
#    sv = SimVis()
#    start_server(sv)
#    asyncio.ensure_future(main(sv))
#    loop = asyncio.get_event_loop()
#    loop.run_forever()

#async def main(sv): #pragma: no cover
#    #create a simple energy
#
#    x = np.linspace(-1, 1, 100)
#    y = x**2
#    sv.set_energy(x, y)
#
#
#
#    while True:
#        sv.set_position(np.random.random(1))
#        await asyncio.sleep(0.5)