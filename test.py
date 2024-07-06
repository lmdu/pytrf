import time
import pytrf
import pyfastx

fa = pyfastx.Fasta('../krait2/data/chr2.fa.gz', uppercase=True)

for s in fa:
	pass

start = time.time()
ssrs = pytrf.GTRFinder(s.name, s.seq).as_list()
print(time.time() - start)
print(len(ssrs))
print(ssrs[0])
print(ssrs[1])

