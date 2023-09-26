def initial_matrix(n, m):
	d = []

	for i in range(n+1):
		d.append([0]*(m+1))
		d[i][0] = i

	for j in range(m+1):
		d[0][j] = j

	return d

def print_matrix(s, ms, mx, st, n, m):
	print('\t\t{}'.format('\t'.join(list(ms))))

	for i in range(n+1):
		if i > 0:
			b = s[st+i]
		else:
			b = ''

		print("{}\t{}".format(b, '\t'.join(map(str, mx[i][0:m+1]))))

def wrap_around_distance(b, s, m, i, d):
	#first pass
	#j = 1
	if b == s[0]:
		c = 0
	else:
		c = 1

	d[i][1] = min(d[i-1][0]+c, d[i-1][m]+c, d[i-1][1]+1)

	#j > 1
	for j in range(2, m+1):
		if b == s[j-1]:
			c = 0
		else:
			c = 1

		d[i][j] = min(d[i-1][j-1]+c, d[i][j-1]+1, d[i-1][j]+1)

	#second pass
	d[i][1] = min(d[i][1], d[i][m]+1)

	for j in range(2, m):
		d[i][j] = min(d[i][j], d[i][j-1]+1)

	return d[i][m] > d[i-1][m]


def wrap_around_extend(s, ms, ml, mx, st, n, me, dr):
	ce = 0

	if not n:
		return 0

	for i in range(1, n+1):
		if wrap_around_distance(s[st+i*dr], ms, ml, i, mx):
			ce += 1

		else:
			ce = 0

		if ce > me:
			break

	i -= ce;

	return i

def wrap_around_backtrace(s, ms, ml, mx, st, i, dr):
	nm = ns = ni = nd = 0
	j = ml

	path = []

	while i > 0 or j > 0:
		path.append((i, j))

		if s[st+i*dr] == ms[j-1]:
			c = 0
		else:
			c = 1

		if j == 1:
			if mx[i][j] == (mx[i-1][ml] + c):
				if c == 0:
					nm += 1
				else:
					ns += 1

				i -= 1
				j = ml

			elif mx[i][j] == (mx[i-1][0] + c):
				if c == 0:
					nm += 1
				else:
					ns += 1

				i -= 1
				j -= 1

			elif mx[i][j] == (mx[i-1][1] + 1):
				ni += 1
				i -= 1

		else:
			"""
			if j < ml and mx[i][j] == (mx[i][j-1] + 1):
				nd += 1
				i -= 1
			
			else:
				if mx[i][j] == (mx[i-1][j-1] + c):
					if c == 0:
						nm += 1
					else:
						ns += 1

					i -= 1
					j -= 1

				elif mx[i][j] == (mx[i][j-1] + 1):
					nd += 1
					j -= 1
				
				elif mx[i][j] == (mx[i-1][j] + 1):
					ni += 1
					i -= 1
			"""

			mval = min(mx[i-1][j-1], mx[i-1][j], mx[i][j-1])

			if mval == mx[i-1][j-1]:
				if mval == mx[i][j]:
					nm += 1
				else:
					ns += 1

				i -= 1
				j -= 1
			elif mval == mx[i-1][j]:
				ni += 1
				i -= 1

			elif mval == mx[i][j-1]:
				nd += 1
				j -= 1

	return nm, ns, ni, nd, path

def atr_finder(seq, max_motif_size=6, min_seed_repeat=3, min_seed_length=10,
				max_consecutive_error=3, min_identity=0.7, max_extend_length=1000):

	matrix = initial_matrix(max_extend_length, max_motif_size)

	size = len(seq)
	atrs = []
	i = 0
	
	while i < size:
		if seq[i] == 'N':
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
					ed = wrap_around_backtrace(seq, motif, j, matrix, extend_start, extend_len, -1)
					tandem_match += ed[0]
					tandem_substitute += ed[1]
					tandem_insert += ed[2]
					tandem_delete += ed[3]

				
				tandem_start = extend_start - extend_len + 1
				
				#extend to right
				extend_start = seed_end
				extend_maxlen = size - extend_start - 1

				if extend_maxlen > max_extend_length:
					extend_maxlen = max_extend_length

				extend_len = wrap_around_extend(seq, motif, j, matrix, extend_start, extend_maxlen,
												max_consecutive_error, 1)

				if extend_len > 0:
					print_matrix(seq, motif, matrix, extend_start, extend_len, j)
					ed = wrap_around_backtrace(seq, motif, j, matrix, extend_start, extend_len, 1)
					tandem_match += ed[0]
					tandem_substitute += ed[1]
					tandem_insert += ed[2]
					tandem_delete += ed[3]

					path = ed[4]

					for a, b in path:
						matrix[a][b] = "{}*".format(matrix[a][b])

					print_matrix(seq, motif, matrix, extend_start, extend_len, j)

					
				
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
	s = "ATGCATGCATGCAGGCTGC"

	atrs = atr_finder(s)

	print(atrs)
