.. note::

	``Stria`` is only responsible for mining tandem repeats from given DNA sequence. It can not directly extract TRs from Fasta files. So we used `pyfastx <https://github.com/lmdu/pyfastx>`_ to parse Fasta file and feed the sequence to ``stria``.

SSR identification
==================

Stria provide ``SSRMiner`` class to find all microsatellites or SSRs from given sequence.

The fastest way to get all SSRs from sequence: 

.. code:: python

	>>> # parse input fasta sequence
	>>> fa = pyfastx.Fastx('tests/data/test.fa.gz', uppercase=True)
	>>> # get the first sequence from fasta
	>>> name, seq, _ = next(fa)

	>>> # feed sequence to stria
	>>> # the fastest way to get all SSRs from sequence
	>>> ssrs = stria.SSRMiner(name, seq).as_list()

You can also iterate over SSRMiner object to get exact tandem repeat (ETR) object. ETR object allows you to access more information and format the information to tsv, csv or gff string.

.. code:: python

	>>> # iterate over ssrminer object to get etr object
	>>> for ssr in stria.SSRMiner(name, seq):
	>>> 	print(ssr.chrom)
	>>> 	print(ssr.motif)
	>>> 	print(ssr.repeats)

You can define the minimum number of repeats required to determine a SSR.

.. code:: python

	>>> # change the minimum repeats for [mono,di,tri,tetra,penta,hexa]
	>>> ssrs = stria.SSRMiner(name, seq, [12,6,4,3,3,3])

A complete example, get all ssrs and output csv format

.. code:: python

	>>> for name, seq, _ in pyfastx.Fastx('tests/data/test.fa', uppercase=True):
	>>> 	for ssr in stria.SSRMiner(name, seq, [12,7,5,4,4,4]):
	>>> 		print(ssr.as_string(','))

VNTR identification
===================

Stria provide ``VNTRMiner`` class to find all minisatellites or VNTRs from given sequence.

The fastest way to get all VNTRs from sequence: 

.. code:: python

	>>> # feed sequence to stria
	>>> # the fastest way to get all VNTRs from sequence
	>>> vntrs = stria.VNTRMiner(name, seq).as_list()

Iterate over vntrminer object to get etr object

.. code:: python

	>>> for vntr in stria.VNTRMiner(name, seq):
	>>> 	print(vntr.chrom)
	>>> 	print(vntr.motif)
	>>> 	print(vntr.repeats)

You can customize the motif size and minimum repeat.

.. code:: python

	>>> vntrs = stria.VNTRMiner(name, seq, min_motif_size=10, max_motif_size=100, min_repeat=3)

A complete example, get all vntrs and output csv format

.. code:: python

	>>> for name, seq, _ in pyfastx.Fastx('tests/data/test.fa', uppercase=True):
	>>> 	for vntr in stria.VNTRMiner(name, seq, 10, 100, 2):
	>>> 		print(vntr.as_string(','))

Exact tandem repeat
===================

When iterating over ``SSRMiner`` or ``VNTRMiner`` object, a exact tandem repeat (ETR) object will be returned.
ETR is a readonly object and allows you to access the attributions and convert to desired formats.

.. code:: python

	>>> ssrs = SSRMiner(name, seq)
	>>> # get one ssr
	>>> ssr = next(ssrs)

	>>> # get sequence name where SSR located on
	>>> ssr.name

	>>> # get one-based start and end position
	>>> ssr.start
	>>> ssr.end

	>>> # get motif sequence
	>>> ssr.motif

	>>> # get number of repeats
	>>> ssr.repeats

	>>> # get length
	>>> ssr.length

	>>> # get SSR sequence
	>>> ssr.seq

	>>> # convert to a list
	>>> ssr.as_list()

	>>> # convert to a dict
	>>> ssr.as_dict()

	>>> # convert to a gff formatted string
	>>> ssr.as_gff()

	>>> # convert to tsv string
	>>> ssr.as_string(separator='\t')

	>>> # convert to csv string
	>>> ssr.as_string(separator=',')

	>>> # added a terminator to the end
	>>> ssr.as_string(separator=',', terminator='\n')

ITR identification
==================

Imperfect tandem repeat
=======================


Commandline interface
=====================