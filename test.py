import time
import pytrf
import pyfastx

fa = pyfastx.Fasta('../data/chr2.fa.gz', uppercase=True)

for s in fa:
	pass

ssrs = pytrf.STRFinder(s.name, s.seq)
for ssr in ssrs:
	print(ssr.as_list())
	break

