#!/usr/bin/env python3

import yt
from yt.funcs import mylog
import numpy as np
import sys
import argparse

parser = argparse.ArgumentParser(
    description='Generate initial conditions for tracer particles.')
parser.add_argument('N', type=int, default=int(1e6),
                    help='Number of particles to draw (default: %(default)s)')
parser.add_argument('-i', '--input', type=str,
                    default='output_00001/info_00001.txt',
                    help='Path to the initial condition '
                    'file (default: %(default)s)')
parser.add_argument('-o', '--output', type=str,
                    default='ic_tracers',
                    help='Name of the output file (ascii, default: '
                    '%(default)s)')
parser.add_argument('--output-type', type=str,
                    default='binary',
                    help='Format of output. Either binary or'
                    'ascii (default: %(default)s)')
parser.add_argument('-x', type=float, nargs=2, default=[0, 1],
                    help='x boundaries')
parser.add_argument('-y', type=float, nargs=2, default=[0, 1],
                    help='y boundaries')
parser.add_argument('-z', type=float, nargs=2, default=[0, 1],
                    help='z boundaries')
parser.add_argument('--filter-by', type=str, default='',
                    help='Filter using this field')
parser.add_argument('--filter-by-min', type=float, default=-np.nan,
                    help='Minimal value (in cgs) for filter')
parser.add_argument('--filter-by-max', type=float, default=np.nan,
                    help='max value (in cgs) for filter')
parser.add_argument('--version', type=int, default=2,
                    help='3 records for version 1 (x, y, z), nparts records for version 2.')
parser.add_argument('--seed', type=int, default=None)

args = parser.parse_args()
if args.seed:
    np.random.seed(args.seed)

yt.enable_parallelism()

Ntracer = args.N
input_fname = args.input
output_fname = args.output

bbox = [args.x[0]]

if args.x[0] > 0 or args.y[0] > 0 or args.y[0] > 0 or \
   args.x[1] < 1 or args.y[1] < 1 or args.z[1] < 1:
    mylog.info('filtering in box')
    left_edge = args.x[0], args.y[0], args.z[0]
    right_edge = args.x[1], args.y[1], args.z[1]
    bbox = np.array([left_edge, right_edge])
    ds = yt.load(input_fname, bbox=bbox)
    ad = ds.box(left_edge=left_edge, right_edge=right_edge)
else:
    mylog.info('getting full box')
    ds = yt.load(input_fname)
    ad = ds.all_data()


# Get position and mass of the cells
x, y, z = xyz = np.array([ad['index', k] for k in 'xyz'])
mass = ad['cell_mass']

# Eventually filter
mask = np.ones(mass.shape, dtype=bool)
if args.filter_by != '':
    val = ad[args.filter_by].in_cgs().value
    if np.isfinite(args.filter_by_min):
        mask &= (val > args.filter_by_min)
    if np.isfinite(args.filter_by_max):
        mask &= (val < args.filter_by_max)

    if mask.sum() == 0:
        mylog.error('No cell found! Aborting')
        print( val.max(), val.min())
        sys.exit(1)

xyzfilter = xyz[:, mask]
massfilter = mass[mask]

# Compute gas mass
gas_mass = massfilter.sum()
mtracer = gas_mass / Ntracer

# Compute number of particle in cell
n_in_cell = massfilter / mtracer

mylog.info(
    'adding particles in cells containing more than 1 particles (%s loops)'
    % np.floor(n_in_cell.max().value).astype(int))
allpos = []
while any(n_in_cell > 1):
    mylog.debug(
        '\tloop: %s' % np.floor(n_in_cell.max().value).astype(int))
    w = (n_in_cell >= 1)

    newpos = xyzfilter[:, w]
    n_in_cell[w] = n_in_cell[w] - 1

    allpos.append(newpos.T)

if len(allpos) > 0:
    pos = np.concatenate(allpos)
else:
    pos = np.empty((0, 3))

mylog.info('adding particles in cells containg more than 0 particles')
d = np.random.rand(n_in_cell.shape[0])
w = d < n_in_cell

newpos = xyzfilter[:, w]
pos = np.concatenate((pos, newpos.T))

mylog.info('saving to %s, %s particles' % (output_fname, pos.shape[0]))
if args.output_type == 'ascii':
    np.savetxt(output_fname, pos, fmt='%.4e %.4e %.4e'.split())
elif args.output_type == 'binary':
    mpart = gas_mass.in_units('code_mass').value / pos.shape[0]

    from scipy.io import FortranFile as FF
    with open(output_fname, 'bw') as f:
        ff = FF(f)
        # Write records
        ff.write_record(pos.shape[0])
        ff.write_record(mpart.astype(np.float64))

        # Version 1: write positions as record
        if args.version == 1:
            ff.write_record(pos[:, 0])
            ff.write_record(pos[:, 1])
            ff.write_record(pos[:, 2])
        else:
            # Version 2: write as binary
            pos.tofile(f)
    print('Set tracer_mass = %s'
          % mpart)
