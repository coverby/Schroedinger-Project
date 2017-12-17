import io
import numpy as np
import Schro1D as schro
import os

def test_read_energy_read():
    '''Tests if the function reads and returns given an input file'''
    test_string = '''
    #test input energy
    #i  V0  c   bs_coefficients   basis_set    domain 
    0,   1,   1,   3,   l,   -inf,    inf
    1,   0,   1,   1,   f,   -inf,    inf
    2,   1,   2,   2,   l,   -1,  1
    3,   1,   2,   2,   f,   -1,  1
    '''


    test_file = 'tests/testfile1.txt'
    indx, V0, const, n_bs, bs, domain = schro.read_energy(test_file)
    print(bs)
    #print(energy)
    #print((np.isclose(pos, [0,1,2,3,4])).any())
    assert((np.isclose(indx, [0,1,2,3])).any())
    assert((np.isclose(V0, [1,0,1,1])).any())
    assert((np.isclose(const,[1,1,2,2])).any())
    assert((np.isclose(n_bs, [3,-1,2,2])).any())
    assert(bs == [b'l', b'f', b'l', b'f'])  #I thought you were my friend, python
    #assert((np.isclose(domain, [-2,-1,4,6,1])).any())