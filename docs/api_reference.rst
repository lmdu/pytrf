API Reference
=============

pytrf.version
-------------

.. py:function:: pytrf.version()

	Get current version of pytrf

	:return: version

	:rtype: str

pytrf.STRFinder
---------------

.. py:class:: pytrf.STRFinder(chrom, seq, mono=12, di=7, tri=5, tetra=4, penta=4, hexa=4)

	Find all exact or perfect short tandem repeats (STRs), simple sequence repeats (SSRs) or microsatellites that meet the minimum repeats on the input sequence

	:param str chrom: the sequence name

	:param str seq: the input DNA sequence

	:param int mono: the minimum tandem repeats for mono-nucleotide repeats

	:param int di: the minimum tandem repeats for di-nucleotide repeats

	:param int tri: the minimum tandem repeats for tri-nucleotide repeats

	:param int tetra: the minimum tandem repeats for tetra-nucleotide repeats

	:param int penta: the minimum tandem repeats for penta-nucleotide repeats

	:param int hexa: the minimum tandem repeats for hexa-nucleotide repeats

	:return: STRFinder object

	.. py:method:: as_list()

		Put all SSRs in a list and return, each SSR in list has 7 columns including [sequence name, start position, end position, motif sequence, motif length, repeats, SSR length]

		:return: all SSRs found

		:rtype: list

pytrf.GTRFinder
---------------

.. py:class:: pytrf.GTRFinder(chrom, seq, max_motif=30, min_repeat=3, min_length=10)

	Find all exact or perfect generic tandem repeats (GTRs) that meet the minimum repeat and minimum length on the input sequence

	:param str chrom: the sequence name

	:param str seq: the input DNA sequence

	:param int max_motif: maximum length of motif sequence

	:param int min_repeat: minimum number of tandem repeats

	:param int min_length: minimum length of tandem repeats

	:return: GTRFinder object

	.. py:method:: as_list()

		Put all GTRs in a list and return, each GTR in list has 7 columns including [sequence name, start position, end position, motif sequence, motif length, repeats, GTR length]

		:return: all GTRs found

		:rtype: list

pytrf.ATRFinder
---------------

.. py:class:: pytrf.ATRFinder(chrom, seq, max_motif_size=6, seed_min_repeat=3, seed_min_length=10, max_continuous_error=3, min_identity=70, max_extend_length=2000)

	Find all approximate or imperfect tandem repeats (ATRs) from the input sequence

	:param str chrom: the sequence name

	:param str seq: the input DNA sequence

	:param int max_motif_size: maximum length of motif

	:param int seed_min_repeat: minimum number of repeat for seed

	:param int seed_min_length: minimum length of seed

	:param int max_continuous_error: maximum number of allowed continuous aligned errors

	:param float min_identity: minimum identity between ATR with its perfect counterpart (0~100)

	:param int max_extend_length: maximum length allowed to extend

	:return: ATRFinder object

	.. py:method:: as_list()

		Put all ATRs in a list and return, each ATR in list has 11 columns including [sequence name, start position, end position, motif sequence, motif length, ATR length, matches, substitutions, insertions, deletions, identity]

pytrf.ETR
---------

.. py:class:: pytrf.ETR

	Readonly exact tandem repeat (ETR) object generated by iterating over STRFinder or GTRFinder object

	.. py:attribute:: chrom

		chromosome or sequence name where ETR located on

	.. py:attribute:: start

		ETR one-based start position on sequence

	.. py:attribute:: end

		ETR one-based end position on sequence

	.. py:attribute:: motif

		motif sequence

	.. py:attribute:: type

		motif length

	.. py:attribute:: repeats

		number of repeats

	.. py:attribute:: length

		length of ETR

	.. py:attribute:: seq

		get the sequence of ETR

	.. py:method:: as_list()

		convert ETR object to a list

	.. py:method:: as_dict()

		convert ETR object to a dict

	.. py:method:: as_gff(terminator='')

		convert ETR object to a gff formatted string

	.. py:method:: as_string(separator='\t', terminator='')

		convert ETR object to a TSV or CSV string by using separator and terminator

		:param str separator: a separator between columns

		:param str terminator: a terminator added to the end of string

		:return: a formatted string

		:rtype: str

pytrf.ATR
---------

.. py:class:: pytrf.ATR

	Readonly imperfect or approximate tandem repeat (ATR) object generated by iterating over ATRFinder object

	.. py:attribute:: chrom

		chromosome or sequence name where ATR located on

	.. py:attribute:: start

		ETR one-based start position on sequence

	.. py:attribute:: end

		ETR one-based end position on sequence

	.. py:attribute:: motif

		motif sequence

	.. py:attribute:: type

		motif length

	.. py:attribute:: length

		length of ITR

	.. py:attribute:: matches

		number of matches

	.. py:attribute:: substitutions

		number of substitutions

	.. py:attribute:: insertions

		number of insertions

	.. py:attribute:: deletions

		number of deletions

	.. py:attribute:: identity

		similar identity

	.. py:attribute:: seq

		get the sequence of ATR

	.. py:method:: as_list()

		convert ATR object to a list

	.. py:method:: as_dict()

		convert ATR object to a dict

	.. py:method:: as_gff(terminator='')

		convert ATR object to a gff formatted string

	.. py:method:: as_string(separator='\t', terminator='')

		convert ATR object to a TSV or CSV string by using separator and terminator

		:param str separator: a separator between columns

		:param str terminator: a terminator added to the end of string

		:return: a formatted string

		:rtype: str