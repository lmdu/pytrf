`Stria` is only responsible for mining tandem repeats from given DNA sequence. It can not directly extract TRs from Fasta files. So we used `pyfastx <https://github.com/lmdu/pyfastx>`_ to parse Fasta file and feed the sequence to `stria`.

SSR identification
==================

Stria provide ``SSRMiner`` class to find all microsatellites from given sequence.

.. cod:: python

	>>> import stria
	>>> import pyfastx

	>>> # parse input fasta sequence
	>>> fa = pyfastx.Fastx('tests/data/test.fa.gz', uppercase=True)
	>>> # get the first sequence from fasta
	>>> name, seq, _ = next(fa)

	>>> # feed sequence to stria
	>>> # the fastest way to get all SSRs from sequence
	>>> ssrs = stria.SSRMiner(name, seq).as_list()

	>>> # iterate over SSRMiner object to get exact tandem repeat (ETR) object
	>>> for ssr in stria.SSRMiner(name, seq):
	>>> 	print(ssr.chrom)
	>>> 	print(ssr.motif)
	>>> 	print(ssr.repeats)

	>>> # change the minimum repeats for [mono,di,tri,tetra,penta,hexa]
	>>> ssrs = stria.SSRMiner(name, seq, [12,6,4,3,3,3])

	>>> # a complete example, get all ssrs and output csv format
	>>> for name, seq, _ in pyfastx.Fastx('tests/data/test.fa', uppercase=True):
	>>> 	for ssr in stria.SSRMiner(name, seq, [12,7,5,4,4,4]):
	>>> 		print(ssr.as_string(','))


VNTR identification
===================

ITR identification
==================


Commandline interface
=====================