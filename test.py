import time
import stripy
import pyfastx

for name, seq, _ in pyfastx.Fastx('tests/data/chr1.fa.gz'):
	start = time.time()
	vntrs = stripy.test(name, seq)
	print(time.time() - start)
