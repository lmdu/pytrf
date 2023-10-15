import pytrf
import pyfastx

fa = pyfastx.Fasta('../krait2/data/chr2.fa.gz', uppercase=True)

for s in fa:
	pass

atrs = pytrf.GTRFinder(s.name, s.seq, min_motif=10, max_motif=100)

for atr in atrs:
	print(atr.as_string())

#print(len(atrs.as_list()))
