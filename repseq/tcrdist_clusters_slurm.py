import argparse
from repseq import tcrdist_clustering as tc


def parse_args():
    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('clonoset_filename', type=str,
                    help='clonoset filename in special tcrdist format')
    parser.add_argument('radius', type=int,
                    help='tcrdist radius')
    parser.add_argument('output_prefix', type=str,
                    help='prefix for output files')
    parser.add_argument('--chain', type=str, choices=['alpha', 'beta'], const="beta",
                    help='chain of TCR: "alpha" or "beta"')
    parser.add_argument('--species', type=str, choices=['human', 'mouse'], const="human",
                    help='species: "human" or "mouse"')
    parser.add_argument('--cpus', type=int, choices=range(1,101), const=1,
                    help='number of cpu cores for parralel calculation')
    parser.add_argument('--group_colname', type=str, const="group",
                    help='number of cpu cores for parralel calculation')

    args = parser.parse_args()
    return args


def main():
    args = parse_args()

    tc.build_tcr_dist_clusters(args.clonoset_filename,
                            args.radius,
                            args.output_prefix,
                            chain=args.chain,
                            species=args.spesies,
                            cpus=args.cpus,
                            group_colname=args.group_colname)


if __name__ == "__main__":
    main()