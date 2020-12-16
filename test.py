import time
import strit
import pyfastx

for name, seq in pyfastx.Fasta('chr1.fa.gz', build_index=False):
	start = time.time()
	ssrs = strit.find_ssrs(seq, (10,7,5,4,4,4))
	print(len(ssrs))
	print(time.time()-start)
