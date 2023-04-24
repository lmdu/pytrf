import time
import stria
import pyfastx
gfile = 'chr1.fa.gz'

for name, seq in pyfastx.Fastx(gfile):
	break

finder = stria.SSRMiner(name, seq)

start = time.time()
ssrs = finder.as_list()
print(time.time() - start)

start = time.time()
ssrs = finder.as_test()
print(time.time() - start)
