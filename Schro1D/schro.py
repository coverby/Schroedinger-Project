import numpy as np
import io
import asyncio
from scipy.integrate import quad

def read_param(input_file):
    #Get the index, potential energy, constant, # of basis set coefficients, basis set, and domain from file
    #Basis set options are Legendre Polynomials (l) or Fourier-Legendre Polynomials (f)
    indx = []
    V0 = []
    const = []
    n_bs = []
    bs = []
    domain = []
    
    data = np.genfromtxt(input_file, dtype=None, delimiter=",",autostrip=True)
    print(data)
    #print(data)
    #print(data[1][1])
    #data2 = pd.read_csv(input_file)
    #print(data2)
    for row in range(len(data)):
        indx.append(data[row][0])
        V0.append(data[row][1])
        const.append(data[row][2])
        n_bs.append(data[row][3])
        bs.append(data[row][4]) 
        domain.append((data[row][5], data[row][6]))
    return indx, V0, const, n_bs, bs, domain


def wavefunc_fou(x, k):
    if (k%2 == 0):
        return np.sin(k*np.pi*x/2)
    else:
        return np.cos(k*np.pi*x/2)

def integrator_fou(fun, k, domain):
    '''Using scipy's quadrature integrator for speed/accuracy/convenience'''
    t0 = domain[0]
    tf = domain[1]
    area = quad(fun, t0, tf, (k))
    return area[0]

def gen_ham(V0, const, n_bs, bs, domain):
    ham = np.zeros([n_bs, n_bs])
    if bs == b'f':
        int_fun = integrator_fou
        fun = wavefunc_fou
        ham[0,0] = 0 + V0*int_fun(fun, 1, domain)  #math index starts at 1, so we adjust k up by 1
        for i in range(1, n_bs):
            ham[i,i] = 1/8*(np.pi**2)*(i**2) + V0*int_fun(fun, i+1, domain)
            for j in range(i-1):
                ham[i,j] = V0*int_fun(fun, j+1, domain)
                ham[j,i] = ham[i,j]
    elif bs == b'l':
        print('Sorry, Legendgre Polynomials not yet supported')
        return 0
    else:
        print('Basis set not recognized!')
        return 0
    return ham


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
    #sv = SimVis()
    print("This is a 1D Numerical Shroedinger's Equation solver")
    parafile = input("Enter the name of the file containing parameter information (format in readme): ")

