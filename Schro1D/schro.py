import numpy as np
import io
import asyncio
from scipy.integrate import quad
from scipy import linalg as sLA
from numpy import linalg as LA

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

def wavefunc_fou2(x,k1,k2):
    if (k1%2 == k2%2):
        if (k1%2 == 0):
            return np.sin(k1*np.pi*x/2)**2
        else:
            return np.cos(k1*np.pi*x/2)**2
    elif (k1%2 == 0):
        return np.sin(k1*np.pi*x/2)*np.cos(k2*np.pi*x/2)
    else:
        return np.cos(k1*np.pi*x/2)*np.sin(k2*np.pi*x/2)

def integrator_fou(fun, k, domain):
    '''Using scipy's quadrature integrator for speed/accuracy/convenience'''
    t0 = domain[0]
    tf = domain[1]
    area = quad(fun, t0, tf, (k))
    return area[0]

def legendre_gen(x, n):
    '''Generate the first n legendre polynomials of x using Bonnet's Recursion'''
    p = []
    if (n < 0):
        return -1
    p.append(1)
    if (n == 1):
        return p
    p.append(x)
    if (n == 2):
        return p
    for i in range(1,n-1):
        p.append(((2*i+1)*x*p[i] - i*p[i-1])/(i+1))
    return p

def legendre_deriv_gen(x, n):
    '''Generates first and second derivatives of the first n legendre polynomials of x'''
    p2 = []
    p2.append(0)
    if (n == 1):
        return p2
    p2.append(0)
    if (n == 2):
        return p2
    p1 = []
    p1.append(0)
    p1.append(1)
    bp = legendre_gen(x,n) #grab the polynomial evaluations at x
    for i in range(2,n):
        p1.append(i*bp[i-1] + x*p1[i-1]) #Very convenient first derivative
        p2.append(i*p1[i-1] + x*p2[i-1] + p1[i-1]) #And easily drived second
    return p2

def legendre_combo(x, n1, n2):
    '''Convenient multiplier for two different legendre polynomials'''
    out1 = legendre_gen(x, n1)[-1]
    out2 = legendre_gen(x, n2)[-1]
    return out1*out2

def legendre_deriv_combo(x, n1, n2):
    '''Combines legendre basis set and legendre second deriv (based on n1) for integration'''
    out1 = legendre_deriv_gen(x, n1)[-1]
    out2 = legendre_gen(x, n2)[-1]
    return out1*out2

def gen_ham(V0, const, n_bs, bs, domain):
    ham = np.zeros([n_bs, n_bs])
    if bs == b'f':
        int_fun = integrator_fou
        fun = wavefunc_fou2
        ham[0,0] = const*1/8*(np.pi*(1))**2 + V0*int_fun(fun, (1,1), domain)  #math index starts at 1?
        for i in range(1, n_bs):
            ham[i,i] = const*1/8*(np.pi*(i+1))**2 + V0*int_fun(fun, (i+1,i+1), domain)
            for j in range(i):
                ham[i,j] = V0*int_fun(fun, (j+1,i+1), domain)
                ham[j,i] = ham[i,j] #orthonormal => symmetric, so we can cheat
    elif bs == b'l':
        int_fun = integrator_fou
        fun1 = legendre_combo
        fun2 = legendre_deriv_combo
        #ham[0,0] = const*int_fun(fun2, (1,1), domain) + V0*int_fun(fun1, (1,1), domain)
        for i in range(0, n_bs):
            ham[i,i] = V0*int_fun(fun1, (i+1,i+1), domain)  #- const*int_fun(fun2, (i+1,i+1), domain)
            for j in range(n_bs):
                ham[i,j] = -const*int_fun(fun2, (i+1,j+1), domain) + V0*int_fun(fun1, (i+1,j+1), domain)
                ham[j,i] = ham[i,j] #orthonormal => symmetric, so we can cheat
    else:
        print('Basis set not recognized!')
        return 0
    return ham

def diagonalize(ham):
    '''Diagonalizes the input hamiltonian and outputs an eigenvector of basis set constants'''
    arr = np.matrix(ham)
    return sLA.eigh(arr) #returns eigenvector w, and (hopefully) ones-diagonal matrix v

def write_output(idx, vout, outname):
    '''Write the output file containing index and basis set coefficients'''
    header_string = "#Index  #Basis Set Coefficients"
    outfile = open(outname, "w")
    outfile.write(header_string)
    outfile.write("\n")
    for i in range(len(idx)):
        outfile.write("\t".join( (str(idx[i]), str(vout[i]) )))
        outfile.write("\n")
    outfile.close()

def main_handler(parafile, outfile): #pragma: no cover
    #This function coordinates the other components after initialization
    indx, V0, const, n_bs, bs, domain = schro.read_param(test_file)
    idxout = []
    vout = []
    for i in range(len(idx)):
        output = diagonalize(gen_ham(V0[i],const[i],n_bs[i],bs[i],domain[i]))

        idxout.append(i)
        vout.append(output)
        
    write_output(idxout,vout,outfile)

def start(): #pragma: no cover
    #sv = SimVis()
    print("This is a 1D Numerical Shroedinger's Equation solver")
    parafile = input("Enter the name of the file containing parameter information (format in readme): ")
    outfile = input("Enter the name of the output file to save output to")
    main_handler()

