import os
import sys
import csv
import stripy
import pyfastx
import argparse
import multiprocessing as mp

def write_to_format(trs, _format, out):
	if _format = 'tsv':
		for tr in trs:
			out.write(tr.as_string('\t', '\n'))

	elif _format = 'csv':
		for tr in trs:
			out.write(tr.as_string(',', '\n'))

	elif _format = 'gff':
		for tr in trs:
			out.write(tr.as_gff())


def find_ssrs(args):
	if args.threads > 1:
		pool = mp.Pool(self.threads)


	else:
		if args.out == 'stdout':
			fw = sys.stdout
		else:
			fw = open(args.out, 'w')

		for name, seq, _ in pyfastx.Fastx(args.fasta):
			ssrs = SSRMiner(name, seq, args.repeats)
			write_to_format(ssrs, args.format, fw)

		if args.out == 'stdout':
			fw.flush()
		else:
			fw.close()


def main():
	parser = argparse.ArgumentParser(
		prog = 'stripy',
		usage = 'stripy COMMAND [OPTIONS]',
		description = "Short tandem repeat sequence identification",
		formatter_class = argparse.RawDescriptionHelpFormatter
	)

	parser.add_argument('-v', '--version',
		action = 'version',
		version = '%(prog)s {}'.format(stripy.version())
	)

	subparsers = parser.add_subparsers(
		title = 'Commands',
		prog = 'stripy',
		metavar = ''
	)

	parser_ssrminer = subparsers.add_parser('ssrminer',
		help = "Find exact microsatellite or simple sequence repeat",
		formatter_class = argparse.ArgumentDefaultsHelpFormatter
	)
	#parser_ssrminer.set_defaults(cmd='')

	parser_ssrminer.add_argument('-r', '--repeats',
		nargs = 6,
		default = [12, 7, 5, 4, 4, 4],
		metavar = ('mono', 'di', 'tri', 'tetra', 'penta', 'hexa'),
		type = int,
		help = "minimum repeats"
	)

	parser_ssrminer.add_argument('-o', '--out',
		default = 'stdout',
		metavar = '',
		type = lambda x: os.path.normcase(x),
		help = "output file"
	)

	parser_ssrminer.add_argument('-f', '--format',
		default = 'tsv',
		choices = ('tsv', 'csv', 'gff'),
		metavar = '',
		help = 'output format, tsv, csv, or gff'
	)

	parser_ssrminer.add_argument('-t', '--threads',
		default = 1,
		type = int,
		metavar = '',
		help = 'number of threads'
	)

	parser_ssrminer.add_argument('fasta',
		help = 'input fasta file, gzip support'
	)


	args = parser.parse_args()



if __name__ == '__main__':
	mp.freeze_support()

	main()
