import argparse


__version__ = '0.1.0'


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Combine chemical shift data.')
    parser.add_argument('outfile', help='Path to output .tsv file.')
    parser.add_argument('--infiles', nargs='+', required=True, help='Paths to .shifts files.')
    parser.add_argument('--efiles', nargs='+', required=True, help='Paths to .energy files.')
    parser.add_argument('--version', '-v', action='version', version=__version__, help='Print version and exit.')

    args = parser.parse_args()

    assert len(args.infiles) == len(args.efiles), 'Number of output and energy files must be equal.'

    args.infiles.sort()
    args.efiles.sort()
