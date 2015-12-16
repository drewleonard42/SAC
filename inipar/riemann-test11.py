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


outname = 'data/riemann-test11_np01.ini'
header0 = 'riemann-test_mhd11\n'
n_dims = 1
dims = ['x']
vars = ['rhop', 'm1', 'Ep', 'b1p', 'rhob', 'Eb', 'b1b']
eqpars = ['gamma', 'eta', 'nu']
eqparvals = [1.4, 0.0, 1.0]

# Define domain
full_domain_size = 256,
full_mincoords = 0.0,
full_maxcoords = 1.0,

# Set uniform values for parameters
varvals = [1.0, 0.0, 2.5, 0.0, 0.0, 0.0, 0.0]

# Get dimensions of process distribution and fiddle some things
procs_index0 = outname.find('_np')+3
nprocs = [int(outname[i:i+2]) for i in range(procs_index0, procs_index0+(2*n_dims), 2)]
domain_size = [int(full_domain_size[i]/nprocs[i]) for i in range(len(full_domain_size))]
mincoords = product(*tuple([np.arange(full_mincoords[dim], full_maxcoords[dim], 1/p) for dim, p in enumerate(nprocs)]))
mincoords = np.array([i for i in mincoords])
maxcoords = mincoords + 1/p

procid = 0
for min, max in zip(mincoords, maxcoords):
    # Define file preamble
    header = header0 + ' {: 6} {: .5E} {: 1} {: 1} {: 1}\n'.format(0, 0.0, n_dims, len(eqpars), len(vars))
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

    # Set initial density on one side of discont.
    outdata[128:, n_dims + vars.index('rhop')] = 0.125

    # Set initial B field on one side of discont.
    outdata[128:, n_dims + vars.index('Ep')] = 0.25
    
    # Output ini info
    ext = '_{:03}.ini'.format(procid)
    np.savetxt(outname.replace('.ini', ext), outdata, header=header, comments="", fmt='% .10E')
    procid += 1
