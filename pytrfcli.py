import sys
import csv
import pytrf
import pyfastx
import argparse
import functools
import multiprocessing

def get_format_result(trs, outfmt, outfw):
	if outfmt == 'tsv':
		for tr in trs:
			print(tr.as_string('\t'), file=outfw)

	elif outfmt == 'csv':
		for tr in trs:
			print(tr.as_string(','), file=outfw)

	elif outfmt == 'gff':
		for tr in trs:
			print(tr.as_gff(), file=outfw)

def str_finder(seq, minrep, outfmt, outfw):
	ssrs = pytrf.STRFinder(seq[0], seq[1], *minrep)
	get_format_result(ssrs, outfmt, outfw)

def gtr_finder(seq, maxmotif, minrep, minlen, outfmt, outfw):
	gtrs = pytrf.GTRFinder(seq[0], seq[1], maxmotif, minrep, minlen)
	get_format_result(gtrs, outfmt, outfw)

def atr_finder(seq, maxmotif, seedrep, seedlen, maxerror, minscore, maxextend, outfmt, outfw):
	atrs = pytrf.ATRFinder(seq[0], seq[1], maxmotif, seedrep, seedlen, maxerror, minscore, maxextend)
	get_format_result(atrs, outfmt, outfw)

def tandem_repeat_finder(args):
	if args.cmd == 'str':
		return functools.partial(str_finder, 
			minrep = args.repeats,
			outfmt = args.out_format,
			outfw = args.out_file
		)

	elif args.cmd == 'gtr':
		return functools.partial(gtr_finder,
			maxmotif = args.max_motif,
			minrep = args.min_repeat,
			minlen = args.min_length,
			outfmt = args.out_format,
			outfw = args.out_file
		)

	elif args.cmd == 'atr':
		return functools.partial(atr_finder,
			maxmotif = args.max_motif_size,
			seedrep = args.min_seed_repeat,
			seedlen = args.min_seed_length,
			maxerror = args.max_continuous_error,
			minscore = args.min_identity,
			maxextend = args.max_extend_length,
			outfmt = args.out_format,
			outfw = args.out_file
		)

def extract_sequence(args):
	fa = pyfastx.Fasta(args.fastx)

	#get input format
	dialect = csv.Sniffer().sniff(args.repeat_file.read(1024))
	args.repeat_file.seek(0)
	reader = csv.reader(args.repeat_file, delimiter=dialect)

	if args.out_format == 'fasta':
		for row in reader:
			chrom = row[0]
			start = int(row[1])
			end = int(row[2])

			left, right = fa.flank(chrom, start, end, args.flank_length, True)
			seq = fa.fetch(chrom, (start, end))
			args.out_file.write(">{}:{}-{}\n{}{}{}".format(chrom, start, end, left, seq, right))

	else:
		if args.out_format == 'csv':
			dialect = ','
		else:
			dialect = '\t'

		writer = csv.writer(args.out_file, delimiter=dialect)

		for row in reader:
			chrom = row[0]
			start = int(row[1])
			end = int(row[2])

			left, right = fa.flank(chrom, start, end, args.flank_length, True)
			seq = fa.fetch(chrom, (start, end))

			row.append(seq)
			row.append(left)
			row.append(right)
			writer.writerow(row)

def main():
	parser = argparse.ArgumentParser(
		prog = 'pytrf',
		usage = 'pytrf command [options] fastx',
		description = "a python package for finding tandem repeats from genomic sequences"
	)

	parser.add_argument('-v', '--version',
		action = 'version',
		version = '%(prog)s {}'.format(pytrf.__version__)
	)

	parser.set_defaults(cmd=None)

	subparsers = parser.add_subparsers(
		title = 'commands',
		prog = 'pytrf',
		metavar = ''
	)

	parser_parent = argparse.ArgumentParser(add_help = False)

	parser_parent.add_argument('fastx',
		metavar = 'fastx',
		help = "input fasta or fastq file (gzip support)"
	)

	parser_parent.add_argument('-o', '--out-file',
		type = argparse.FileType('wt'),
		default = sys.stdout,
		metavar = '',
		help = "output file (default: stdout)"
	)

	parser_parent.add_argument('-f', '--out-format',
		default = 'tsv',
		metavar = '',
		help = "output format, tsv, csv or gff (default: tsv)"
	)

	#ssr finder
	parser_ssrfinder = subparsers.add_parser('findstr',
		help = "find exact or perfect short tandem repeats",
		parents = [parser_parent]
	)
	parser_ssrfinder.set_defaults(cmd='str')

	parser_ssrfinder.add_argument('-r', '--repeats',
		nargs = 6,
		default = [12, 7, 5, 4, 4, 4],
		metavar = ('mono', 'di', 'tri', 'tetra', 'penta', 'hexa'),
		type = int,
		help = "minimum repeats for each STR type (default: 12 7 5 4 4 4)"
	)

	#gtr finder
	parser_gtrfinder = subparsers.add_parser('findgtr',
		help = "find exact or perfect generic tandem repeats",
		parents = [parser_parent]
	)
	parser_gtrfinder.set_defaults(cmd='gtr')

	parser_gtrfinder.add_argument('-m', '--max-motif',
		default = 30,
		metavar = '',
		type = int,
		help = "maximum motif length (default: 30)"
	)

	parser_gtrfinder.add_argument('-r', '--min-repeat',
		default = 3,
		metavar = '',
		type = int,
		help = "minimum repeat number (default: 3)"
	)

	parser_gtrfinder.add_argument('-l', '--min-length',
		default = 10,
		metavar = '',
		type = int,
		help = "minimum repeat length (default: 10)"
	)

	#atr finder
	parser_atrfinder = subparsers.add_parser('findatr',
		help = "find approximate or imperfect tandem repeats",
		parents = [parser_parent]
	)
	parser_atrfinder.set_defaults(cmd='atr')

	parser_atrfinder.add_argument('-m', '--max-motif-size',
		default = 6,
		metavar = '',
		type = int,
		help = "maximum motif length (default: 6)"
	)

	parser_atrfinder.add_argument('-r', '--min-seed-repeat',
		default = 3,
		metavar = '',
		type = int,
		help = "minimum repeat number for seed (default: 3)"
	)

	parser_atrfinder.add_argument('-l', '--min-seed-length',
		default = 10,
		metavar = '',
		type = int,
		help = "minimum length for seed (default: 10)"
	)

	parser_atrfinder.add_argument('-e', '--max-continuous-error',
		default = 3,
		metavar = '',
		type = int,
		help = "maximum number of continuous alignment errors (default: 3)"
	)

	parser_atrfinder.add_argument('-p', '--min-identity',
		default = 70,
		metavar = '',
		type = float,
		help = "minimum identity from 0 to 100 (default: 70)"
	)

	parser_atrfinder.add_argument('-x', '--max-extend-length',
		default = 2000,
		metavar = '',
		type = int,
		help = "maximum length allowed to extend (default: 2000)"
	)

	#extract sequence
	parser_extract = subparsers.add_parser('extract',
		help = "get tandem repeat sequence and flanking sequence",
		#conflict_handler = 'resolve',
		#parents = [parser_parent]
	)
	parser_extract.set_defaults(cmd='get')

	parser_extract.add_argument('fastx',
		metavar = 'fastx',
		help = "input fasta or fastq file (gzip support)"
	)

	parser_extract.add_argument('-r', '--repeat-file',
		metavar = '',
		required = True,
		type = argparse.FileType('rt'),
		help = "the csv or tsv output file of findatr, findstr or findgtr"
	)

	parser_extract.add_argument('-o', '--out-file',
		type = argparse.FileType('wt'),
		default = sys.stdout,
		metavar = '',
		help = "output file (default: stdout)"
	)

	parser_extract.add_argument('-f', '--out-format',
		default = 'tsv',
		metavar = '',
		help = "output format, tsv, csv or fasta (default: tsv)"
	)

	parser_extract.add_argument('-l', '--flank-length',
		default = 100,
		metavar = '',
		type = int,
		help = "flanking sequence length (default: 100)"
	)

	args = parser.parse_args()

	if args.cmd is None:
		parser.print_help()

	elif args.cmd == 'get':
		with args.out_file, args.repeat_file:
			extract_sequence(args)

	else:
		with args.out_file:
			func = tandem_repeat_finder(args)
			seqs = pyfastx.Fastx(args.fastx)

			for seq in seqs:
				func(seq)
