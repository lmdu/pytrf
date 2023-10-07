import pytrf
import pyfastx

fa = pyfastx.Fasta('../krait2/data/chr2.fa.gz', uppercase=True)

for s in fa:
	pass

atrs = pytrf.ATRFinder(s.name, s.seq)

for atr in atrs:
	print(atr.as_string())
