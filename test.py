import pytrf
import pyfastx

for s in pyfastx.Fasta('../krait2/data/chr2.fa.gz'):
	pass

finder = pytrf.ATRFinder(s.name, s.seq[0:13965])

n = 0

for atr in finder:
	print(atr.as_list())

print('####################')

atrs = finder.as_list()

for atr in atrs:
	print(atr)

