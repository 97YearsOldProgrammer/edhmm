import isoform2
import argparse

parser = argparse.ArgumentParser(
	description='test gtag sites')
parser.add_argument('fasta')
parser.add_argument('--min_intron', required=False, type=int, default=35,
	metavar='<int>', help='minimum length of intron [%(default)i]')
parser.add_argument('--min_exon', required=False, type=int, default=25,
	metavar='<int>', help='minimum length exon [%(default)i]')
parser.add_argument('--flank', required=False, type=int, default=99,
	metavar='<int>', help='genomic flank on each side [%(default)i]')
arg = parser.parse_args()

name, seq  = next( isoform2.read_fasta(arg.fasta) )
dons, accs = isoform2.gtag_sites(seq, arg.flank, arg.min_exon)

print(dons)
print(accs)