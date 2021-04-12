import time
import stripy
import pyfastx

for name, seq, _ in pyfastx.Fastx('tests/data/chr1.fa.gz'):
	vntrs = stripy.VNTRMiner(name, seq, 7, 100, 3)
	for vntr in vntrs:
		print(vntr.as_string())
