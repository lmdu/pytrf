import time
import stripy
import pyfastx

for name, seq in pyfastx.Fasta('tests/data/chr1.fa.gz', build_index=False, uppercase=True):
	vntrs = stripy.ITRMiner(name, seq).as_list()
	for vntr in vntrs:
		print("\t".join(map(str, vntr)))
