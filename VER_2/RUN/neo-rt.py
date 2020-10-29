#!/usr/bin/python3
from glob import glob
import sys
sys.path.append('../../../libneo/python')
from ninjarun import read_config, build_and_run

target = 'neo-rt'
sources = [
    'main.f90',
    'common.f90',
    'orbit_{orbit_type}.f90',
    'transport_{transport_type}.f90'
]

config = read_config('neo-rt.in')

if config['orbit_type'] == 'full':
    sources += glob('../../POTATO/VER_1/SRC/*.f')
    sources += glob('../../POTATO/VER_1/SRC/*.f90')
    sources.remove('../../POTATO/VER_1/SRC/tt.f90')

build_and_run(target, sources, config=config, fflags='-lblas -llapack')
