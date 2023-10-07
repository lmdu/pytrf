import sys

def initial_matrix(n, m):
	d = []

	for i in range(n+1):
		d.append([0]*(m+1))
		d[i][0] = i

	for j in range(m+1):
		d[0][j] = j

	return d

def print_matrix(seq, mseq, matrix, start, n, mlen, direction=1):
	print('\t\t{}'.format('\t'.join(list(mseq))))

	for i in range(n+1):
		if i > 0:
			base = seq[start+i*direction]
		else:
			base = ''

		print("{}\t{}".format(base, '\t'.join(map(str, matrix[i][0:mlen+1]))))

def wrap_around_distance(base, mseq, mlen, i, matrix):
	#first pass
	#j = 1
	if base == mseq[0]:
		cost = 0
	else:
		cost = 1

	matrix[i][1] = min(matrix[i-1][0]+cost, matrix[i-1][mlen]+cost, matrix[i-1][1]+1)

	#j > 1
	for j in range(2, mlen+1):
		if base == mseq[j-1]:
			cost = 0
		else:
			cost = 1

		matrix[i][j] = min(matrix[i-1][j-1]+cost, matrix[i][j-1]+1, matrix[i-1][j]+1)

	#second pass
	#j = 1
	matrix[i][1] = min(matrix[i][1], matrix[i][mlen]+1)

	#j > 1
	for j in range(2, mlen):
		matrix[i][j] = min(matrix[i][j], matrix[i][j-1]+1)

	return matrix[i][mlen] > matrix[i-1][mlen]


def wrap_around_extend(seq, mseq, mlen, matrix, start, size, max_error, direction):
	current_error = 0

	if size <= 0: return 0

	for i in range(1, size+1):
		if wrap_around_distance(seq[start+i*direction], mseq, mlen, i, matrix):
			current_error += 1

		else:
			current_error = 0

		if current_error > max_error: break

	i -= current_error

	return i

def wrap_around_backtrace(mlen, matrix, i):
	num_mat = num_sub = num_ins = num_del = 0
	j = mlen

	path = []

	while i > 0 or j > 0:
		path.append((i, j))

		if j == 1:
			v = min(matrix[i][mlen], matrix[i-1][mlen], matrix[i-1][0], matrix[i-1][1])

			if i > 0 and j < mlen and v == matrix[i][mlen]:
				num_del += 1
				j = mlen

			elif v == matrix[i-1][mlen]:
				if v == matrix[i][j]:
					num_mat += 1
				else:
					num_sub += 1

				i -= 1
				j = mlen

			elif v == matrix[i-1][0]:
				if v == matrix[i][j]:
					num_mat += 1
				else:
					num_sub += 1

				i -= 1
				j -= 1

			elif v == matrix[i-1][1]:
				num_ins += 1
				i -= 1

		else:
			v = min(matrix[i-1][j-1], matrix[i-1][j], matrix[i][j-1])

			if v == matrix[i-1][j-1]:
				if v == matrix[i][j]:
					num_mat += 1
				else:
					num_sub += 1

				i -= 1
				j -= 1
			elif v == matrix[i-1][j]:
				num_ins += 1
				i -= 1

			elif v == matrix[i][j-1]:
				num_del += 1
				j -= 1

	return num_mat, num_sub, num_ins, num_del, path

def atr_finder(seq, max_motif_size=6, min_seed_repeat=3, min_seed_length=10,
				max_consecutive_error=3, min_identity=0.7, max_extend_length=1000):

	matrix = initial_matrix(max_extend_length, max_motif_size)

	size = len(seq)
	atrs = []
	i = 0

	while i < size:
		if seq[i] == 'N':
			i += 1
			continue

		seed_start = i

		for j in range(1, max_motif_size+1):
			b = size - j;

			while i < b and seq[i] == seq[i+j]:
				i += 1

			seed_length = i + j - seed_start
			seed_repeat = int(seed_length / j)

			seed_length = seed_repeat * j

			if seed_repeat >= min_seed_repeat and seed_length >= min_seed_length:
				motif = seq[seed_start:seed_start+j]

				seed_end = seed_start + seed_length - 1
				tandem_match = seed_length
				tandem_substitute = 0
				tandem_insert = 0
				tandem_delete = 0

				#extend to left
				extend_start = seed_start
				extend_maxlen = extend_start

				if extend_maxlen > max_extend_length:
					extend_maxlen = max_extend_length

				extend_len = wrap_around_extend(seq, motif[::-1], j, matrix, extend_start,
												 extend_maxlen, max_consecutive_error, -1)

				if extend_len > 0:
					print("left:")
					print_matrix(seq, motif, matrix, extend_start, extend_len, j, -1)
					ed = wrap_around_backtrace(j, matrix, extend_len)
					tandem_match += ed[0]
					tandem_substitute += ed[1]
					tandem_insert += ed[2]
					tandem_delete += ed[3]

					#path = ed[4]

					#for a, b in path:
					#	matrix[a][b] = "{}*".format(matrix[a][b])
				
				tandem_start = extend_start - extend_len + 1
				
				#extend to right
				extend_start = seed_end
				extend_maxlen = size - extend_start - 1

				if extend_maxlen > max_extend_length:
					extend_maxlen = max_extend_length

				extend_len = wrap_around_extend(seq, motif, j, matrix, extend_start, extend_maxlen,
												max_consecutive_error, 1)

				if extend_len > 0:
					#print_matrix(seq, motif, matrix, extend_start, extend_len, j)
					#ed = wrap_around_backtrace(j, matrix, extend_len)
					tandem_match += ed[0]
					tandem_substitute += ed[1]
					tandem_insert += ed[2]
					tandem_delete += ed[3]

					path = ed[4]

					#for a, b in path:
					#	matrix[a][b] = "{}*".format(matrix[a][b])

					#print_matrix(seq, motif, matrix, extend_start, extend_len, j)
				
				tandem_align = tandem_match + tandem_insert + tandem_substitute + tandem_delete
				tandem_identity = tandem_match / tandem_align

				if tandem_identity >= min_identity:
					tandem_end = extend_start + extend_len + 1
					tandem_length = tandem_end - tandem_start + 1

					atrs.append((motif, j, tandem_start, tandem_end, tandem_length, tandem_match,
								tandem_substitute, tandem_insert, tandem_delete, tandem_identity))

					i = tandem_end
					break

			i = seed_start

		i += 1

	return atrs

if __name__ == '__main__':
	#s = "AAGAAGAAGAAGCCGAGAAGGTAGATAG"
	#s = "ATGCATGCATGCAGGCTGC"

	import pyfastx
	for s in pyfastx.Fasta('../krait2/data/chr2.fa.gz'):
		pass
	atrs = atr_finder(s.seq)
