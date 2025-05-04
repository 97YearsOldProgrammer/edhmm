import os
import glob
import argparse

import isoform2
import hints

parser = argparse.ArgumentParser(
	description='Alternative isoform generator')
parser.add_argument('dir', type=str, metavar='<fastas dir>',
	help='dir for the small gene set')
parser.add_argument('geniso', type=str, metavar='<geniso exe>',
    help='input geniso ext path')
parser.add_argument('model', type=str, metavar='<splice model>',
	help='input splicemodel file')
parser.add_argument('hmm', type=str, metavar='<hmm path>',
	help='input hmm exe path')
parser.add_argument('--min_intron', required=False, type=int, default=35,
	metavar='<int>', help='minimum length of intron [%(default)i]')
parser.add_argument('--min_exon', required=False, type=int, default=25,
	metavar='<int>', help='minimum length exon [%(default)i]')
parser.add_argument('--flank', required=False, type=int, default=99,
	metavar='<int>', help='genomic flank on each side [%(default)i]')
parser.add_argument('--limit_fasta', required=False, type=int, default=1,
	metavar='<int>', help='limit number of fasta [%(default)i]')
parser.add_argument('--limit_isoform', required=False, type=int, default=10000,
	metavar='<int>', help='boundary for isoform [%(default)i]')
arg = parser.parse_args()

def run(dir):
    fastas = glob.glob(os.path.join(dir, "*.fa"))
    fcount = 0
     
    for fasta in fastas:
        if fcount > arg.limit_fasta: break
         
        with open(fasta, 'r') as file:
            # get donors and acceptors first
            name, seq  = next(isoform2.read_fasta(fasta))
            dons, accs = isoform2.gtag_sites(seq, arg.flank, arg.min_exon)
            
            # isoform limit check
            count = hints.countiso(dons, accs, arg.min_intron, arg.min_exon, arg.limit_isoform)
            if count < arg.limit_isoform: continue
            fcount += 1
                         
            # run geniso
            output  = hints.run_geniso2(arg.geniso, fasta, arg.model) 
            output2 = hints.run_geniso2(arg.geniso, fasta, arg.model, arg.hmm) 

run(arg.dir)
