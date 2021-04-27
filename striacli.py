import os
import sys
import csv
import time
import queue
import stria
import shutil
import pyfastx
import argparse
import multiprocessing as mp

def format_and_write_to_file(args, out, tres):
	if args.out_format == 'tsv':
		for tre in tres:
			out.write(tre.as_string('\t', '\n'))

	elif args.out_format == 'csv':
		for tre in tres:
			out.write(tre.as_string(',', '\n'))

	elif args.out_format == 'gff':
		for tre in tres:
			out.write(tre.as_gff())

def find_tandem_repeats_by_type(name, seq, args):
	if args.cmd == 'ssr':
		tres = stria.SSRMiner(name, seq, args.repeats)

	elif args.cmd == 'vntr':
		tres = stria.VNTRMiner(name, seq, args.min_motif_size, 
								args.max_motif_size, args.min_repeats)
	elif args.cmd == 'itr':
		tres = stria.ITRMiner(name, seq, args.min_motif_size, args.max_motif_size,
								args.seed_min_repeats, args.seed_min_length,
								args.max_continuous_errors, args.substitution_penalty,
								args.insertion_penalty, args.deletion_penalty,
								args.min_match_ratio, args.max_extend_size)

	return tres


def find_tandem_repeats_worker(args, tasks, event, worker_id):
	with open("{}.{}".format(args.out_file, worker_id), 'w') as fw:
		while 1:
			if event.is_set() and tasks.empty():
				break

			try:
				name, seq = tasks.get_nowait()
			except:
				time.sleep(0.01)
				continue

			tres = find_tandem_repeats_by_type(name, seq, args)
			format_and_write_to_file(args, fw, tres)

def find_tandem_repeats_with_multicore(args):
	pool = mp.Pool(args.threads)
	manager = mp.Manager()
	tasks = manager.Queue(args.threads*2)
	event = manager.Event()

	#add workers
	l = len(str(args.threads))
	for i in range(args.threads):
		work_id = str(i).zfill(l)
		pool.apply_async(find_tandem_repeats_worker, (args, tasks, event, work_id))

	#add tasks
	for name, seq, _ in pyfastx.Fastx(args.fasta, uppercase=True):
		tasks.put((name, seq), block=True, timeout=None)

	event.set()
	pool.close()
	pool.join()

	#merge results
	with open(args.out_file, 'wb') as fw:
		for i in range(args.threads):
			temp_file = "{}.{}".format(args.out_file, str(i).zfill(l))

			if os.path.isfile(temp_file):
				with open(temp_file, 'rb') as fh:
					shutil.copyfileobj(fh, fw)

				#remove temp file
				os.remove(temp_file)

def find_tandem_repeats_with_singlecore(args):
	with open(args.out_file, 'w') as fw:
		for name, seq, _ in pyfastx.Fastx(args.fasta, uppercase=True):
			tres = find_tandem_repeats_by_type(name, seq, args)
			format_and_write_to_file(args, fw, tres)

def find_tandem_repeats(args):
	if args.threads > 1:
		find_tandem_repeats_with_multicore(args)
	else:
		find_tandem_repeats_with_singlecore(args)

def main():
	parser = argparse.ArgumentParser(
		prog = 'stria',
		usage = 'stria COMMAND [OPTIONS]',
		description = "short tandem repeat identification and analysis",
		formatter_class = argparse.RawDescriptionHelpFormatter
	)

	parser.add_argument('-v', '--version',
		action = 'version',
		version = '%(prog)s {}'.format(stria.version())
	)

	subparsers = parser.add_subparsers(
		title = 'Commands',
		prog = 'stria',
		metavar = ''
	)

	#ssr miner
	parser_ssrminer = subparsers.add_parser('ssrminer',
		help = "Find exact microsatellites or simple sequence repeats",
		formatter_class = argparse.ArgumentDefaultsHelpFormatter
	)
	parser_ssrminer.set_defaults(cmd='ssr')
	parser_ssrminer.set_defaults(func=find_tandem_repeats)

	parser_ssrminer.add_argument('-r', '--repeats',
		nargs = 6,
		default = [12, 7, 5, 4, 4, 4],
		metavar = ('mono', 'di', 'tri', 'tetra', 'penta', 'hexa'),
		type = int,
		help = "minimum repeats"
	)

	parser_ssrminer.add_argument('-o', '--out-file',
		default = 'stria_ssrs.out',
		metavar = '',
		type = lambda x: os.path.normcase(x),
		help = "output file"
	)

	parser_ssrminer.add_argument('-f', '--out-format',
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

	#vntr miner
	parser_vntrminer = subparsers.add_parser('vntrminer',
		help = "Find exact minisatellites or variable number tandem repeats",
		formatter_class = argparse.ArgumentDefaultsHelpFormatter
	)
	parser_vntrminer.set_defaults(cmd='vntr')
	parser_vntrminer.set_defaults(func=find_tandem_repeats)

	parser_vntrminer.add_argument('-m', '--min-motif-size',
		default = 7,
		metavar = '',
		type = int,
		help = "minimum motif length"
	)

	parser_vntrminer.add_argument('-M', '--max-motif-size',
		default = 30,
		metavar = '',
		type = int,
		help = "maximum motif length"
	)

	parser_vntrminer.add_argument('-r', '--min-repeats',
		default = 2,
		metavar = '',
		type = int,
		help = "minimum repeat number"
	)

	parser_vntrminer.add_argument('-o', '--out-file',
		default = 'stria_vntrs.out',
		metavar = '',
		type = lambda x: os.path.normcase(x),
		help = "output file"
	)

	parser_vntrminer.add_argument('-f', '--out-format',
		default = 'tsv',
		choices = ('tsv', 'csv', 'gff'),
		metavar = '',
		help = 'output format, tsv, csv, or gff'
	)

	parser_vntrminer.add_argument('-t', '--threads',
		default = 1,
		type = int,
		metavar = '',
		help = 'number of threads'
	)

	parser_vntrminer.add_argument('fasta',
		help = 'input fasta file, gzip support'
	)

	#itr miner
	parser_itrminer = subparsers.add_parser('itrminer',
		help = "Find imperfect tandem repeats",
		formatter_class = argparse.ArgumentDefaultsHelpFormatter
	)
	parser_itrminer.set_defaults(cmd='itr')
	parser_itrminer.set_defaults(func=find_tandem_repeats)

	parser_itrminer.add_argument('-m', '--min-motif-size',
		default = 1,
		metavar = '',
		type = int,
		help = "minimum motif length"
	)

	parser_itrminer.add_argument('-M', '--max-motif-size',
		default = 6,
		metavar = '',
		type = int,
		help = "maximum motif length"
	)

	parser_itrminer.add_argument('-r', '--seed-min-repeats',
		default = 3,
		metavar = '',
		type = int,
		help = "minimum repeat number for seed"
	)

	parser_itrminer.add_argument('-l', '--seed-min-length',
		default = 8,
		metavar = '',
		type = int,
		help = "minimum length for seed"
	)

	parser_itrminer.add_argument('-e', '--max-continuous-errors',
		default = 2,
		metavar = '',
		type = int,
		help = "maximum number of continuous alignment errors"
	)

	parser_itrminer.add_argument('-s', '--substitution-penalty',
		default = 0.5,
		metavar = '',
		type = float,
		help = "substitution penalty"
	)

	parser_itrminer.add_argument('-i', '--insertion-penalty',
		default = 1.0,
		metavar = '',
		type = float,
		help = "insertion penalty"
	)

	parser_itrminer.add_argument('-d', '--deletion-penalty',
		default = 1.0,
		metavar = '',
		type = float,
		help = "deletion penalty"
	)

	parser_itrminer.add_argument('-p', '--min-match-ratio',
		default = 0.7,
		metavar = '',
		type = float,
		help = "extending match ratio"
	)

	parser_itrminer.add_argument('-x', '--max-extend-size',
		default = 2000,
		metavar = '',
		type = int,
		help = "maximum length allowed to extend"
	)

	parser_itrminer.add_argument('-o', '--out-file',
		default = 'stria_itrs.out',
		metavar = '',
		type = lambda x: os.path.normcase(x),
		help = "output file"
	)

	parser_itrminer.add_argument('-f', '--out-format',
		default = 'tsv',
		choices = ('tsv', 'csv', 'gff'),
		metavar = '',
		help = 'output format, tsv, csv, or gff'
	)

	parser_itrminer.add_argument('-t', '--threads',
		default = 1,
		type = int,
		metavar = '',
		help = 'number of threads'
	)

	parser_itrminer.add_argument('fasta',
		help = 'input fasta file, gzip support'
	)

	args = parser.parse_args()
	args.func(args)


if __name__ == '__main__':
	#mp.freeze_support()
	main()
