Pytrf
#####

.. image:: https://github.com/lmdu/pytrf/actions/workflows/wheel.yml/badge.svg
   :target: https://github.com/lmdu/pytrf/actions/workflows/wheel.yml
   :alt: Github Action

.. image:: https://readthedocs.org/projects/pytrf/badge/?version=latest
   :target: https://pytrf.readthedocs.io/en/latest/?badge=latest
   :alt: Readthedocs

.. image:: https://img.shields.io/pypi/v/pytrf.svg
   :target: https://pypi.org/project/pytrf
   :alt: PyPI

.. image:: https://img.shields.io/pypi/pyversions/pytrf
	:target: https://pypi.org/project/pytrf
	:alt: PyPI

.. image:: https://app.codacy.com/project/badge/Grade/bbe59e55f686465ca5824c69583e9718
	:target: https://app.codacy.com/gh/lmdu/pytrf/dashboard?utm_source=gh&utm_medium=referral&utm_content=&utm_campaign=Badge_grade

*a fast Python package for finding tandem repeat sequences*

Introduction
============

A Tandem repeat (TR) in genomic sequence is a set of adjacent short DNA
sequence repeated consecutively. The pytrf is a lightweight Python C extension for identification of tandem repeats. The pytrf enables to fastly identify both exact
or perfect SSRs. It also can find generic tandem repeats with any size of motif,
such as with maximum motif length of 100 bp. Additionally, it has capability of finding approximate or imperfect tandem repeats. Furthermore, the pytrf not only can be used as Python package but also provides command line interface for users to facilitate the identification of tandem repeats.

Usage
=====

The pytrf can be used as Python package. It requires `pyfastx <https://github.com/lmdu/pyfastx>`_ to parse FASTA or FASTQ file.

.. code:: python

	>>> import pytrf
	>>> import pyfastx
	>>> fa = pyfastx.Fastx('test.fa', uppercase=True):
	>>> for name, seq in fa:
	>>> 	for ssr in STRFinder(name, seq):
	>>> 		print(ssr.as_string())

Command line
============

The pytrf also provides command line tools for you to find tandem repeats from given FASTA or FASTQ file.

.. code:: sh

	pytrf -h

	usage: pytrf command [options] fastx

	a python package for finding tandem repeats from genomic sequences

	options:
	  -h, --help     show this help message and exit
	  -v, --version  show program's version number and exit

	commands:

	    findstr      find exact or perfect short tandem repeats
	    findgtr      find exact or perfect generic tandem repeats
	    findatr      find approximate or imperfect tandem repeats
	    extract      get tandem repeat sequence and flanking sequence

For example:

.. code:: sh

	pytrf findstr test.fa

Documentation
=============

For more detailed usage, see our manual: `https://pytrf.readthedocs.io <https://pytrf.readthedocs.io>`_
