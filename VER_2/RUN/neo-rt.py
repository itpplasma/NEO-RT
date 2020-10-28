#!/usr/bin/python
# see https://gist.github.com/syoyo/778d2294fd5534e0f923
from os import getcwd, chdir, path, popen, system
import sys
import f90nml

fc = 'gfortran'
target = 'neo-rt.x'
sources = [
    'main.f90', 'orbit_{orbit_type}.f90', 'transport_{transport_type}.f90'
]

origdir = getcwd()
rundir = path.dirname(path.abspath(__file__))
builddir = path.abspath(path.join(rundir, '..', 'BUILD'))
srcdir = path.abspath(path.join(rundir, '..', 'SRC'))
sys.path.append(srcdir)
import ninja_syntax

config = f90nml.read('neo-rt.in')['config']
sourcepaths = []
for k in range(len(sources)):
    sourcepaths.append(path.join('..', 'SRC', sources[k].format_map(config)))

status_build = -1
status_run = -1
try:
    chdir(builddir)
    with open('build.ninja', 'w') as ninjafile:
        ninja = ninja_syntax.Writer(ninjafile)
        ninja.rule('fc', f'{fc} $in -o $out')
        ninja.build(target, 'fc', sourcepaths)
    status_build = system('ninja')

finally:
    chdir(origdir)

if status_build == 0:  # build has succeeded
    try:
        #status_run = system(path.join(builddir, target))
        pipe = popen(path.join(builddir, target))
        print(pipe.read(), end='')
        status_run = pipe.close()
    finally:
        if status_run is not None:
            print(f'Error: {target} exited with code {status_run}.')
else:
    print(f'Error: Exiting due to build error for {target}.')
