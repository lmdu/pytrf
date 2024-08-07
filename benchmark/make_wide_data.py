import sys

print("tool\tgenome\tmemory\ttime")
with open(sys.argv[1]) as fh:
	header = fh.readline().split('\t')[3:]
	for line in fh:
		cols = line.strip().split('\t')

		for idx, col in enumerate(cols[3:]):
			if idx % 2 == 0:
				tool = header[idx]
				memory = round(float(col)/1024, 2)
			else:
				if tool not in ['PERF', 'RPTRF', 'Mreps']:
					print(tool, cols[0], memory, col, sep='\t')
