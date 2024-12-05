Welcome to pytrf's documentation!
==================================

A Tandem repeat (TR) in genomic sequence is a set of adjacent short DNA
sequence repeated consecutively. The core sequence or repeat unit is generally
called motif. According to the motif length, tandem repeats can be classified
as microsatellites and minisatellites. Microsatellites are also known as simple
sequence repeats (SSRs) or short tandem repeats (STRs) with motif length of 1-6 bp.
Minisatellites are also sometimes referred to as variable number of tandem repeats
(VNTRs) has longer motif length than micorsatellites.

The pytrf is a lightweight Python C extension for identification of tandem repeats.
The pytrf enables to fastly identify both exact or perfect SSRs. It also can find generic
tandem repeats with any size of motif, such as with maximum motif length of 100 bp.
Additionally, it has capability of finding approximate or imperfect tandem repeats.
Furthermore, the pytrf not only can be used as Python package but also provides command
line interface for users to facilitate the identification of tandem repeats.

Note: pytrf is not a Python binding to common used tool `TRF <https://tandem.bu.edu/trf/trf.html>`_.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   installation
   usage
   changelog
   api_reference


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
