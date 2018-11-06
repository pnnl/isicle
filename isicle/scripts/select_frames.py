import argparse
from os.path import *
import pandas as pd
from isicle import __version__


def select_frames(path, nframes=10, low=1.25E6, high=1.45E6):
    tmp = []
    with open(path, 'r') as f:
        for line in f:
            if 'NSTEP' in line:
                props = [x for x in line.split() if x is not '=']
                d = {'step': int(props[1]),
                     'time': float(props[3]),
                     'temp': float(props[5])}
                tmp.append(d)

    df = pd.DataFrame(tmp)
    stepsize = df['step'][1] - df['step'][0]
    df['frame'] = df['step'] // stepsize

    ss = df[(df['step'] >= low) & (df['step'] <= high)].sample(n=nframes)

    return ss['frame'].values


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Select frames from simulated annealing')
    parser.add_argument('infile', help='path to simulated annealing .out file')
    parser.add_argument('crdfile', help='path to .crd file')
    parser.add_argument('outfiles', nargs='+', help='paths to output .trajin files (one per frame)')
    parser.add_argument('--nframes', type=int, default=10, help='number of frames')
    parser.add_argument('--low', type=float, default=1.25e+06, help='lower timestep bound')
    parser.add_argument('--high', type=float, default=1.45e+06, help='upper timestep bound')
    parser.add_argument('-v', '--version', action='version', version=__version__, help='print version and exit')

    args = parser.parse_args()

    frames = select_frames(args.infile,
                           nframes=args.nframes,
                           low=args.low,
                           high=args.high)

    for frame, outfile in zip(frames, args.outfiles):
        with open(outfile, 'w') as f:
            f.write('trajin %s %s %s\n' % (args.crdfile, frame, frame))
            f.write('trajout %s.mol2 mol2\n' % splitext(outfile)[0])
