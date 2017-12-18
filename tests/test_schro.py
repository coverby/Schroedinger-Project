import io
import numpy as np
import Schro1D as schro
import os
import math

def test_read_param_read():
    '''Tests if the function reads and returns given an input file'''
    test_string = '''
    #test paramter read
    #i  V0  c   bs_coefficients   basis_set    domain 
    0,   1,   1,   3,   l,   -inf,    inf
    1,   0,   1,   1,   f,   -inf,    inf
    2,   1,   2,   2,   l,   -1,  1
    3,   1,   2,   2,   f,   -1,  1
    '''


    test_file = 'tests/testfile1.txt'
    indx, V0, const, n_bs, bs, domain = schro.read_param(test_file)
    #print((np.isclose(pos, [0,1,2,3,4])).any())
    assert((np.isclose(indx, [0,1,2,3])).any())
    assert((np.isclose(V0, [1,0,1,1])).any())
    assert((np.isclose(const,[1,1,2,2])).any())
    assert((np.isclose(n_bs, [3,-1,2,20])).any())
    assert(bs == [b'l', b'f', b'l', b'f'])  #I thought you were my friend, python
    assert((np.isclose(domain, [(-math.inf, math.inf), (-math.inf, math.inf), (-1, 1), (-1, 1)])).any())

def test_wavefunc_fou():
    '''Quick sanity check for the Fourier basis wavefunction'''
    assert(np.isclose(schro.wavefunc_fou(0,1), 1))
    assert(np.isclose(schro.wavefunc_fou(0,0), 0))

def test_math_integrator():
    '''Performs a quick test on the integrator function on sin(x)'''
    tstart = 0
    tstop = 2
    k = 1
    yout = schro.integrator_fou(schro.wavefunc_fou,k, [tstart, tstop])    
    assert(np.isclose(yout,0))

def test_ham_generator():
    '''Tests the hamiltonian generator.  Unfortunately, products are not edible'''
    V0 = 1
    const = 1
    n_bs = 5
    bs = b'f'
    domain = [-1,1]
    ham = np.matrix(schro.gen_ham(V0, const, n_bs, bs, domain))
    print(ham)
    #The Hamiltonian should be self-adjunct so...
    assert(np.isclose(ham.H, ham).any())