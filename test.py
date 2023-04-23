import sys
import time
import pyfastx

start_time = time.time()
mins = [0, 12, 14, 15, 16, 20, 24]

for name, seq in pyfastx.Fastx(sys.argv[1]):
	pass

size = len(seq)
count = 0
i = 0
jump = 0

while i < size:
	if seq[i] == 'N':
		i += 1
		continue

	start = i
	run = [0, 0, 0, 0, 0, 0, 0]

	for j in range(1, 7):
		'''
		if j > 1:
			if run[1] >= j:
				jump += 1
				continue

			if j == 4 and run[2] >= j:
				jump += 1
				continue

			elif j == 6 and (run[2] >= j or run[3] >= j):
				jump += 1
				continue
		'''

		b = size - j

		while i < b and seq[i] == seq[i+j]:
			i += 1

		run[j] = i + j - start

		if run[j] >= mins[j]:
			count += 1
			repeat = run[j]//j
			i = start + repeat*j
			run = [0, 0, 0, 0, 0, 0, 0]
			break
		else:
			i = start

	i += 1

end_time = time.time()

print(end_time-start_time)
print(count)
print(jump)
print(size)
