#! /usr/bin/env python2
from __future__ import print_function

import sys
sys.path.append("../../../")

from fabricate import *

try:
    from scripts import sacconfig
    cfg = sacconfig.SACConfig()
except ImportError:
    cfg = None

if cfg:
    F_compiler = cfg.compiler
    F_flags = cfg.compiler_flags.split(' ')
    options = cfg.vac_modules
else:
    print("No config module found using default values for mpi on intel")
    F_compiler = 'mpif90'
    F_flags = ['-free', '-mcmodel=medium']
    options = ['vaccd', 'vacmpi']

F_ext = '.f90'
pre_processor = './vacpp.sh'
sources = ['vac', 'vacio', 'vacgrid', 'vacphys0', 'vacphys', 'vacusr']
includes = [ 'vacdef', 'vacpar', 'vacusrpar']

def build():
    pre_process()
    compile_modules()
    compile()
    link()

def pre_process():
    for source in sources + includes + options:
        run(pre_processor, source+'.t',  source+F_ext)

def compile_modules():
    for source in includes:
        run(F_compiler, *F_flags + ['-c', source+F_ext])

def compile():
    for source in sources + options:
        run(F_compiler, *F_flags + ['-c', source+F_ext])

def link():
    objects = [s+'.o' for s in sources + options]
    run(F_compiler, '-o', '../vac', objects)

def clean():
    autoclean()

if __name__ == '__main__':
    main()
