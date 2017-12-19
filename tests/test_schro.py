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

def test_wavefunc_fou2():
    '''Quick sanity check for the two-k Fourier basis wavefunction'''
    assert(np.isclose(schro.wavefunc_fou2(0,1,0), 0))
    assert(np.isclose(schro.wavefunc_fou2(0,0,0), 0))
    assert(np.isclose(schro.wavefunc_fou2(0,1,1), 1))

def test_legendre_gen():
    '''Test the output of the recursive legendre polynomial generator'''
    n = 4
    x1 = 0
    x2 = 1
    x3 = -1
    print(schro.legendre_gen(x1, n))
    print(schro.legendre_gen(x2, n))
    assert(np.isclose(schro.legendre_gen(x1, n), [1, 0, -.5, 0]).all())
    assert(np.isclose(schro.legendre_gen(x2, n), [1, 1, 1, 1]).all())
    assert(np.isclose(schro.legendre_gen(x3, n), [1, -1, 1, -1]).all())
    #This is such a nice basis set for these sorts of things

def test_legendre_deriv_gen():
    '''Test the output of the legendre polynomial second derivative generator'''
    n = 5
    x1 = 0
    x2 = 1
    x3 = -1
    print(schro.legendre_deriv_gen(x1, n))
    print(schro.legendre_deriv_gen(x2, n))
    print(schro.legendre_deriv_gen(x3, n))
    assert(np.isclose(schro.legendre_deriv_gen(x1, n), [0, 0, 3, 0, -7.5]).all())
    assert(np.isclose(schro.legendre_deriv_gen(x2, n), [0, 0, 3, 15, 45]).all())
    assert(np.isclose(schro.legendre_deriv_gen(x3, n), [0, 0, 3, -15, 45]).all())
    #This is such a nice basis set for these sorts of things

def test_math_integrator():
    '''Performs a quick test on the integrator function on sin(x)'''
    tstart = 0
    tstop = 2
    k = 7
    yout = schro.integrator_fou(schro.wavefunc_fou,k, [tstart, tstop])  
    print(yout)  
    assert(np.isclose(yout,0))

def test_math_integrator2():
    '''Performs a quick test on the integrator function on fourier basis'''
    tstart = -1
    tstop = 1
    k1 = 3
    k2 = 2
    yout = schro.integrator_fou(schro.wavefunc_fou2,(k1, k2), [tstart, tstop])  
    print(yout)  
    assert(np.isclose(yout,0))
    k1 = 2
    k2 = 2
    yout = schro.integrator_fou(schro.wavefunc_fou2,(k1, k2), [tstart, tstop])  
    print(yout)  
    assert(np.isclose(yout,1))

def test_math_integrator_leg():
    '''Performs a quick test on the integrator function legendre polynomials'''
    tstart = -1
    tstop = 1
    n1  = 2
    n2 = 3
    yout = schro.integrator_fou(schro.legendre_combo, (n1, n2), [tstart, tstop])  
    print(yout)  
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

def test_diagonalize():
    '''Tests the hamiltonian diagonalizer via the normalized eigenvectors '''
    V0 = 1
    const = 1
    n_bs = 6
    bs = b'f'
    domain = [-1,1]
    ham = np.matrix(schro.gen_ham(V0, const, n_bs, bs, domain))
    #The Hamiltonian should be self-adjunct so...
    eig, v = schro.diagonalize(ham)
    #print(ham)
    #print(eig)
    #print(v)
    #print(np.diagonal(np.absolute(v)).sum())
    assert(np.isclose(np.diagonal(np.absolute(v)).sum(), n_bs,rtol=1.e-2))