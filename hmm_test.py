import os
import glob
import argparse

import isoform2
import hints

parser = argparse.ArgumentParser(
	description='Alternative isoform generator')
parser.add_argument('dir', type=str, metavar='<fastas dir>',
	help='dir for the small gene set')
parser.add_argument('hmm', type=str, metavar='<hmm path>',
	help='input hmm exe path')
parser.add_argument('--min_intron', required=False, type=int, default=35,
	metavar='<int>', help='minimum length of intron [%(default)i]')
parser.add_argument('--min_exon', required=False, type=int, default=25,
	metavar='<int>', help='minimum length exon [%(default)i]')
parser.add_argument('--flank', required=False, type=int, default=99,
	metavar='<int>', help='genomic flank on each side [%(default)i]')
arg = parser.parse_args()

def run(dir):
    fastas = glob.glob(os.path.join(dir, "*.fa"))
    for fasta in fastas:         
        with open(fasta, 'r') as file:
            # test whether the ouput from model is correct
            name, seq    = next(isoform2.read_fasta(fasta))
            dons1, accs1 = hints.mgtag_sites(seq, arg.flank, arg.min_exon, 39)
            dons2, accs2 = isoform2.hmm_sites(arg.hmm, fasta)
            dons = [don[0] for don in dons2]
            accs = [acc[0] for acc in accs2]
            print(dons1)
            print(dons)
            print(accs)
            print(accs1)
            assert all([
    		    len(dons1) == len(dons),
    		    len(accs1) == len(accs),
    		    dons1 == dons,
    		    accs1 == accs
		    ]), "Check hmm output. Splice site not match"

run(arg.dir)