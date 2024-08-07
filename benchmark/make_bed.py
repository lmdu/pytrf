import os
import sys

REPS = [0, 13, 7, 5, 4, 3, 3]

def make_bed_for_gmata(infile, outdir):
	outfile = os.path.join(outdir, os.path.basename(infile).replace('.out', '.bed'))
	with open(infile) as fh, open(outfile, 'w') as fw:
		fh.readline()
		for line in fh:
			cols = line.strip().split()

			chrom = cols[0].strip('>').split()[0]
			if int(cols[4]) >= REPS[len(cols[5])]:
				print(chrom, cols[2], cols[3], sep='\t', file=fw)

def make_bed_for_kmerssr(infile, outdir):
	outfile = os.path.join(outdir, os.path.basename(infile).replace('.out', '.bed'))
	with open(infile) as fh, open(outfile, 'w') as fw:
		fh.readline()
		for line in fh:
			cols = line.strip().split()

			chrom = cols[0].split()[0]
			start = int(cols[3]) + 1
			end = int(cols[3]) + len(cols[1]) * int(cols[2])

			print(chrom, start, end, sep='\t', file=fw)

def make_bed_for_misa(infile, outdir):
	outfile = os.path.join(outdir, os.path.basename(infile).replace('.out', '.bed'))
	with open(infile) as fh, open(outfile, 'w') as fw:
		fh.readline()

		for line in fh:
			cols = line.strip().split()
			chrom = cols[0].split()[0]
			print(chrom, cols[5], cols[6], sep='\t', file=fw)

def make_bed_for_phobos(infile, outdir):
	outfile = os.path.join(outdir, os.path.basename(infile).replace('.out', '.bed'))
	with open(infile) as fh, open(outfile, 'w') as fw:
		for line in fh:
			if line[0] == '#':
				continue

			if line.startswith('seq-name'):
				continue

			cols = line.strip().split()
			chrom = cols[0].split()[0]
			if  float(cols[6]) >= REPS[int(cols[1])]:
				print(chrom, cols[3], cols[4], sep='\t', file=fw)

def make_bed_for_sciroko(infile, outdir):
	outfile = os.path.join(outdir, os.path.basename(infile).replace('.td', '.bed'))
	with open(infile) as fh, open(outfile, 'w') as fw:
		for line in fh:
			if line.startswith('Seq_Name'):
				continue

			if line[0] == '#':
				continue

			if not line.strip():
				continue

			cols = line.strip().split()
			chrom = cols[0].split()[0]
			print(chrom, cols[3], cols[4], sep='\t', file=fw)

def make_bed_for_ssrit(infile, outdir):
	outfile = os.path.join(outdir, os.path.basename(infile).replace('.out', '.bed'))
	with open(infile) as fh, open(outfile, 'w') as fw:
		for line in fh:
			cols = line.strip().split()
			chrom = cols[0].split()[0]
			print(chrom, cols[5], cols[6], sep='\t', file=fw)

def make_bed_for_tantan(infile, outdir):
	outfile = os.path.join(outdir, os.path.basename(infile).replace('.out', '.bed'))
	with open(infile) as fh, open(outfile, 'w') as fw:
		for line in fh:
			cols = line.strip().split()

			if float(cols[4]) >= REPS[int(cols[3])]:
				chrom = cols[0].split()[0]
				print(chrom, cols[1], cols[2], sep='\t', file=fw)

def make_bed_for_other(infile, outdir):
	outfile = os.path.join(outdir, os.path.basename(infile).replace('.out', '.bed'))
	with open(infile) as fh, open(outfile, 'w') as fw:
		for line in fh:
			cols = line.strip().split()

			print(cols[0], cols[1], cols[2], 1, sep='\t', file=fw)

if __name__ == '__main__':
	outdir = sys.argv[1]
	infiles = sys.argv[2:]

	for infile in infiles:
		fname = os.path.basename(infile)
		tool = fname.split('.')[1]

		if tool == 'gmata':
			make_bed_for_gmata(infile, outdir)
		elif tool == 'kmer':
			make_bed_for_kmerssr(infile, outdir)
		elif tool == 'misa':
			make_bed_for_misa(infile, outdir)
		elif tool == 'phobos':
			make_bed_for_phobos(infile, outdir)
		elif tool == 'sciroko':
			make_bed_for_sciroko(infile, outdir)
		elif tool == 'ssrit':
			make_bed_for_ssrit(infile, outdir)
		elif tool == 'tantan':
			make_bed_for_tantan(infile, outdir)
		else:
			make_bed_for_other(infile, outdir)
