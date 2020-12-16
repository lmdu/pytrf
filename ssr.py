import pyfastx
minreps = [0, 12, 7, 5, 4, 4, 4]

for _, seq, _ in pyfastx.Fastx('chr1.fa.gz'):
	pass

size = len(seq)
lens = []
ssrs = []

i = 0
while (i < size):
	if seq[i] == 'N':
		continue

	j = 1
	while j <= 6:
		start = i
		while (i+j) < size and seq[i] == seq[i+j]:
			i += 1

		replen = i + j - start
		repeats = replen/j;

		if repeats >= minreps[j]:
			motif = seq[start:start+j]
			i += j
			ssrs.append((motif, j, repeats, start+1, i, replen))

			j = 0

		else:
			i = start

		j += 1

	i += 1
