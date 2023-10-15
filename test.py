import pytrf
import pyfastx

fa = pyfastx.Fasta('../krait2/data/chr2.fa.gz', uppercase=True)

for s in fa:
	pass

atrs = pytrf.ATRFinder(s.name, s.seq, min_motif_size=10, max_motif_size=100, min_seed_repeat=2, max_consecutive_error=3)

for atr in atrs:
	atr.as_string()

#print(len(atrs.as_list()))
