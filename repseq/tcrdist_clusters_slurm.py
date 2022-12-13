import argparse

import sys
import os
REPSEQ_PATH = os.path.join(os.path.expanduser("~"), "soft", "repseq")
sys.path.append(REPSEQ_PATH)
from repseq import tcrdist_clustering as tc



def parse_args():
    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('clonoset_filename', type=str,
                    help='clonoset filename in special tcrdist format')
    parser.add_argument('radius', type=int,
                    help='tcrdist radius')
    parser.add_argument('output_prefix', type=str,
                    help='prefix for output files')
    parser.add_argument('--chain', type=str, choices=['alpha', 'beta'], default="beta",
                    help='chain of TCR: "alpha" or "beta"')
    parser.add_argument('--species', type=str, choices=['human', 'mouse'], default="human",
                    help='species: "human" or "mouse"')
    parser.add_argument('--cpus', type=int, choices=range(1,101), default=1,
                    help='number of cpu cores for parralel calculation')
    parser.add_argument('--group_colname', type=str, default="group",
                    help='number of cpu cores for parralel calculation')
    parser.add_argument('--cdr3_gap_penalty', type=int, default=4,
                    help='cdr3 region gap penalty')

    args = parser.parse_args()
    return args


def main():
    
    args = parse_args()
    print("Read args: {} {} {} {} {} {} {} {}".format(args.clonoset_filename,args.radius,
                                                   args.output_prefix,args.chain,
                                                   args.species,args.cpus,args.group_colname,args.cdr3_gap_penalty))
    tc.build_tcr_dist_clusters(args.clonoset_filename,
                            args.radius,
                            args.output_prefix,
                            chain=args.chain,
                            species=args.species,
                            cpus=args.cpus,
                            group_colname=args.group_colname,
                            cdr3_gap_penalty=args.cdr3_gap_penalty)


if __name__ == "__main__":
    main()