Usage
=====

.. note::

	``Pytrf`` is only responsible for finding tandem repeats from given DNA sequence. It can not directly extract TRs from Fasta files. So we used `pyfastx <https://github.com/lmdu/pyfastx>`_ to parse Fasta file and feed the sequence to ``pytrf``.

STR identification
------------------

Pytrf provide ``STRFinder`` class to find all microsatellites or SSRs from given sequence.

The fastest way to get all SSRs from sequence: 

.. code:: python

	>>> # parse input fasta sequence
	>>> fa = pyfastx.Fastx('tests/data/test.fa.gz', uppercase=True)
	>>> # get the first sequence from fasta
	>>> name, seq = next(fa)

	>>> # feed sequence to pytrf
	>>> # the fastest way to get all SSRs from sequence
	>>> ssrs = pytrf.STRFinder(name, seq).as_list()

You can also iterate over STRFinder object to get exact tandem repeat (ETR) object. ETR object allows you to access more information and format the information to tsv, csv or gff string.

.. code:: python

	>>> # iterate over STRFinder object to get ETR object
	>>> for ssr in pytrf.STRFinder(name, seq):
	>>> 	print(ssr.chrom)
	>>> 	print(ssr.motif)
	>>> 	print(ssr.repeat)

You can define the minimum number of repeats required to determine a SSR.

.. code:: python

	>>> # change the minimum repeats for mono-, di-, tri-, tetra-, penta-, hexa-nucleotide repeat
	>>> ssrs = pytrf.STRFinder(name, seq, 10, 6, 4, 3, 3, 3)

A complete example, get all ssrs and output csv format

.. code:: python

	>>> fa = pyfastx.Fastx('tests/data/test.fa', uppercase=True)
	>>> for name, seq in fa:
	>>> 	for ssr in pytrf.STRFinder(name, seq, 12, 7, 5, 4, 4, 4):
	>>> 		print(ssr.as_string(','))

GTR identification
-------------------

Pytrf provide ``GTRFinder`` class to find all generic tandem repeats (GTRs) with
any size of motif from given sequence.

The fastest way to get all GTRs from sequence:

.. code:: python

	>>> # feed sequence to pytrf
	>>> # the fastest way to get all GTRs from sequence
	>>> gtrs = pytrf.GTRFinder(name, seq).as_list()

Iterate over GTRFinder object to get ETR object

.. code:: python

	>>> for gtr in pytrf.GTRMiner(name, seq):
	>>> 	print(gtr.chrom)
	>>> 	print(gtr.motif)
	>>> 	print(gtr.repeat)

You can customize the motif size, minimum repeat and minimum length.

.. code:: python

	>>> gtrs = pytrf.GTRFinder(name, seq, min_motif=20, max_motif=100, min_repeat=3, min_length=10)

A complete example, get all gtrs and output csv format

.. code:: python

	>>> fa = pyfastx.Fastx('tests/data/test.fa', uppercase=True):
	>>> for name, seq in fa:
	>>> 	for gtr in pytrf.GTRFinder(name, seq, 30, 100, 2, 10):
	>>> 		print(vntr.as_string(','))

Exact tandem repeat
-------------------

When iterating over ``STRFinder`` or ``GTRFinder`` object, an exact tandem repeat (ETR) object will be returned.
ETR is a readonly object and allows you to access the attributes and convert to desired formats.

.. code:: python

	>>> ssrs = STRFinder(name, seq)
	>>> # get one ssr
	>>> ssr = next(ssrs)

	>>> # get sequence name where SSR located on
	>>> ssr.chrom

	>>> # get one-based start and end position
	>>> ssr.start
	>>> ssr.end

	>>> # get repeat sequence
	>>> ssr.seq

	>>> # get motif sequence
	>>> ssr.motif

	>>> # get number of repeats
	>>> ssr.repeat

	>>> # get repeat length
	>>> ssr.length

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

ATR identification
------------------

Pytrf provide ``ATRFinder`` class to find all imperfect or approximate tandem repeats from given sequence.

The fastest way to get all ATRs from sequence:

.. code:: python

	>>> # feed sequence to pytrf
	>>> # the fastest way to get all ATRs from sequence
	>>> itrs = pytrf.ATRFinder(name, seq).as_list()

Iterate over ATRFinder object to get atr object

.. code:: python

	>>> for atr in pytrf.ATRFinder(name, seq):
	>>> 	print(atr.chrom)
	>>> 	print(atr.motif)
	>>> 	print(atr.length)

You can customize the motif size and seed parameters.

.. code:: python

	>>> itrs = pytrf.ATRFinder(name, seq, max_motif_size=10, seed_min_repeat=3, seed_min_length=10)

A complete example, get all atrs and output csv format

.. code:: python

	>>> fa = pyfastx.Fastx('tests/data/test.fa', uppercase=True)
	>>> for name, seq in fa:
	>>> 	for atr in pytrf.ATRFinder(name, seq):
	>>> 		print(atr.as_string(','))

Approximate tandem repeat
-------------------------

When iterating over ``ATRFinder`` object, an imperfect or approximate tandem repeat (ATR) object will be returned.
ATR is a readonly object and allows you to access the attributes and convert to desired formats.

.. code:: python

	>>> atrs = ATRFinder(name, seq)
	>>> # get one ATR
	>>> atr = next(atrs)

	>>> # get sequence name where ATR located on
	>>> atr.name

	>>> # get one-based start and end position
	>>> atr.start
	>>> atr.end

	>>> # get repeat sequence
	>>> atr.seq

	>>> # get motif sequence
	>>> atr.motif

	>>> # get length
	>>> atr.length

	>>> # get number of matches
	>>> atr.matches

	>>> # get number of substitutions
	>>> atr.substitutions

	>>> # get number of insertions
	>>> atr.insertions

	>>> # get number of deletions
	>>> atr.deletions

	>>> # convert to a list
	>>> atr.as_list()

	>>> # convert to a dict
	>>> atr.as_dict()

	>>> # convert to a gff formatted string
	>>> atr.as_gff()

	>>> # convert to tsv string
	>>> atr.as_string(separator='\t')

	>>> # convert to csv string
	>>> atr.as_string(separator=',')

	>>> # added a terminator to the end
	>>> atr.as_string(separator=',', terminator='\n')

Commandline interface
---------------------

``Pytrf`` also provide command line tools for users to find tandem repeats from fasta or fastq files.

.. code:: sh

	pytrf -h

	usage: pytrf command [options] fastx

	a python package for finding tandem repeats from genomic sequences

	options:
	  -h, --help     show this help message and exit
	  -v, --version  show program version number and exit

	commands:

	    findstr      find exact or perfect short tandem repeats
	    findgtr      find exact or perfect generic tandem repeats
	    findatr      find approximate or imperfect tandem repeats
	    extract      get tandem repeat sequence and flanking sequence

Find exact microsatellites or simple sequence repeats (SSRs) from fasta/q file.

.. code:: sh

	pytrf findstr -h

	usage: pytrf findstr [-h] [-o] [-f] [-r mono di tri tetra penta hexa] fastx

	positional arguments:
	  fastx                 input fasta or fastq file (gzip support)

	options:
	  -h, --help            show this help message and exit
	  -o , --out-file       output file (default: stdout)
	  -f , --out-format     output format, tsv, csv or gff (default: tsv)
	  -r mono di tri tetra penta hexa, --repeats mono di tri tetra penta hexa
	                        minimum repeats for each STR type (default: 12 7 5 4 4 4)

Find exact generic tandem repeats (GTRs) from fasta/q file.

.. code:: sh

	pytrf gtrfinder -h

	usage: pytrf findgtr [-h] [-o] [-f] [-m] [-M] [-r] [-l] fastx

	positional arguments:
	  fastx               input fasta or fastq file (gzip support)

	options:
	  -h, --help          show this help message and exit
	  -o , --out-file     output file (default: stdout)
	  -f , --out-format   output format, tsv, csv or gff (default: tsv)
	  -m , --min-motif    minimum motif length (default: 10)
	  -M , --max-motif    maximum motif length (default: 100)
	  -r , --min-repeat   minimum repeat number (default: 3)
	  -l , --min-length   minimum repeat length (default: 10)

Find imperfect or approximate tandem repeats (ATRs)

.. code:: sh

	pytrf atrfinder -h

	usage: pytrf findatr [-h] [-o] [-f] [-m] [-M] [-r] [-l] [-e] [-p] [-x] fastx

	positional arguments:
	  fastx                 input fasta or fastq file (gzip support)

	options:
	  -h, --help            show this help message and exit
	  -o , --out-file       output file (default: stdout)
	  -f , --out-format     output format, tsv, csv or gff (default: tsv)
	  -m , --min-motif-size
	                        minimum motif length (default: 1)
	  -M , --max-motif-size
	                        maximum motif length (default: 6)
	  -r , --min-seed-repeat
	                        minimum repeat number for seed (default: 3)
	  -l , --min-seed-length
	                        minimum length for seed (default: 10)
	  -e , --max-continuous-error
	                        maximum number of continuous alignment errors (default: 3)
	  -p , --min-identity   minimum identity from 0 to 1 (default: 0.7)
	  -x , --max-extend-length
	                        maximum length allowed to extend (default: 2000)

Extract tandem repeat sequence and flanking sequence according results of findatr, findgtr or findstr.

.. code:: sh

	pytrf extract -h

	usage: pytrf extract [-h] -r  [-o] [-f] [-l] fastx

	positional arguments:
	  fastx                 input fasta or fastq file (gzip support)

	options:
	  -h, --help            show this help message and exit
	  -r , --repeat-file    the csv or tsv output file of findatr, findstr or findgtr
	  -o , --out-file       output file (default: stdout)
	  -f , --out-format     output format, tsv, csv or fasta (default: tsv)
	  -l , --flank-length   flanking sequence length (default: 100)
