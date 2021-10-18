import time
import stria
import pyfastx

start = time.time()
for name, seq, _ in pyfastx.Fastx('hg19.fa'):
	ssrs = stria.test(name, seq)
	#ssrs = stria.SSRMiner(name,seq).as_list()
print(time.time()-start)
#print(len(ssrs))
