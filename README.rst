stria
#####

.. image:: https://github.com/lmdu/stria/actions/workflows/main.yml/badge.svg
   :target: https://github.com/lmdu/stria/actions/workflows/main.yml
   :alt: Github Action

.. image:: https://readthedocs.org/projects/stria/badge/?version=latest
   :target: https://stria.readthedocs.io/en/latest/?badge=latest
   :alt: Readthedocs

.. image:: https://img.shields.io/pypi/v/stria.svg
   :target: https://pypi.org/project/stria
   :alt: PyPI

*a fast Python package for identification and analysis of short tandem repeat sequences*

Introduction
============

A Tandem repeat (TR) in genomic sequence is a set of adjacent short DNA sequence repeated consecutively.
The core sequence or repeat unit is generally called motif. According to the motif length, tandem repeats
can be classified as microsatellites and minisatellites. Microsatellites are also known as simple sequence
repeats (SSRs) or short tandem repeats (STRs) with motif length of 1-6 bp. Minisatellites are also sometimes
referred to as variable number of tandem repeats (VNTRs) has longer motif length than micorsatellites.

The ``stria`` is a lightweight Python C extension for identification and analysis of short tandem repeats.
The ``stria`` enables to fastly identify both exact and imperfect SSRs and VNTRs from large numbers of DNA sequences.
The ``stria`` also provides command line tools for users to extract tandem repeats from Fasta files.

Usage
=====

.. code:: python

	>>> import stria
	>>> import pyfastx
	>>> for name, seq, _ in pyfastx.Fastx('test.fa.gz'):
	>>> 	for ssr in SSRMiner(name, seq):
	>>> 		print(ssr.as_string())

Documentation
=============

`https://stria.readthedocs.io/ <https://stria.readthedocs.io/>`_
