from __future__ import division
import numpy as np
from itertools import product
import sys


def prod(sequence):
    """Quick and dirty function to return integer size of domain for given dimensions"""
    product = 1
    for x in sequence:
        product *= x
    return product


outname = 'data/orszagtang-test22_np0202.ini'
header0 = 'orszagtang-test_mhd22\n'
n_dims = 2
n_dims2 = 2
dims = ['x', 'y']
vars = ['rhop', 'm1', 'm2', 'Ep', 'b1p', 'b2p', 'rhob', 'Eb', 'b1b', 'b2b']
eqpars = ['gamma', 'eta']
eqparvals = [5.0/3.0, 0.0]

# Define domain
full_domain_size = 256, 256
full_mincoords = 0.0, 0.0
full_maxcoords = 1.0, 1.0

# Set uniform values for parameters
varvals = [25/(36*np.pi), 0.0, 0.0, 0.0, 0.0, 0.0]

# Get dimensions of process distribution and fiddle some things
procs_index0 = outname.find('_np')+3
nprocs = [int(outname[i:i+2]) for i in range(procs_index0, procs_index0+(2*n_dims), 2)]
domain_size = [int(full_domain_size[i]/nprocs[i]) for i in range(len(full_domain_size))]
mincoords = product(*tuple([np.arange(full_mincoords[dim], full_maxcoords[dim], 1/p) for dim, p in enumerate(nprocs)])) # This will break if {x|y|z}_max != 1
mincoords = np.array([i for i in mincoords])
maxcoords = mincoords + 1/p # As above

procid = 0
for min, max in zip(mincoords, maxcoords):
    # Define file preamble
    header = header0 + ' {: 6} {: .5E} {: 1} {: 1} {: 1}\n'.format(0, 0.0, n_dims, n_dims2, len(vars))
    for x in domain_size: header += ' {}'.format(x)
    header += '\n'
    for x in eqparvals: header += ' {: .5E}'.format(x)
    header += '\n'
    for x in dims + vars + eqpars: header += '{} '.format(x)

    # Arrange values in array for output
    outdata = np.zeros(shape=(prod(domain_size), len(dims + vars)))
    coords = product(*tuple([np.arange(min[i], max[i], (max[i]-min[i])/domain_size[i]) for i in range(len(dims))]))
    coords = np.array([i for i in coords])
    for i, dim in enumerate(dims):
        coords[:, i] += (0.5 * ((max[i]-min[i]) / domain_size[i]))
    ##### This next bit may or may not be a great big hack #####
    coords = coords[:, ::-1]
    ##### <\hack>
    outdata[:, :n_dims] = coords
    for i, val in enumerate(varvals):
        outdata[:, n_dims+i] = val
    """
     Manually adjust values as needed
    """
    # Set initial B field
    B0 = 1 / np.sqrt(4*np.pi)
    outdata[:, n_dims + vars.index('b1p')] = -B0*np.sin(4*np.pi*outdata[:, dims.index('y')])
    outdata[:, n_dims + vars.index('b2p')] = B0*np.sin(2*np.pi*outdata[:, dims.index('x')])
    
    # Set initial velocity
    rho = varvals[vars.index('rhop')]
    outdata[:, n_dims + vars.index('m1')] = np.sin(2*np.pi*outdata[:, dims.index('y')]) * rho
    outdata[:, n_dims + vars.index('m2')] = np.sin(2*np.pi*outdata[:, dims.index('x')]) * rho
    
    # Set initial thermal energy
    P = 5/(12*np.pi)
    gamma = eqparvals[eqpars.index('gamma')]
    vsq = ((outdata[:, n_dims + vars.index('m1')] / rho)**2) + ((outdata[:, n_dims + vars.index('m2')] / rho)**2)
    Bsq = (outdata[:, n_dims + vars.index('b1p')]**2) + (outdata[:, n_dims + vars.index('b2p')]**2)
    outdata[:, n_dims + vars.index('Ep')] = (P/(gamma-1)) + (0.5*(rho*vsq)) + (0.5*Bsq)
    
    # Output ini info
    ext = '_{:03}.ini'.format(procid)
    np.savetxt(outname.replace('.ini', ext), outdata, header=header, comments="")#, fmt='% .10E')
    procid += 1
