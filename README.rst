minerva_barcoded_read_deconvolution
===================================

.. image:: https://img.shields.io/pypi/v/minerva_barcoded_read_deconvolution.svg
    :target: https://pypi.python.org/pypi/minerva_barcoded_read_deconvolution
    :alt: Latest PyPI version


This is a demonstration program for graph-based deconvolution of linked reads. A more efficient multi-threaded version of this software is currently being developed. This program was used to write 'Minerva: An Alignment and Reference Free Approach to Deconvole Linked-Reads for Metagenomics'


Usage
-----

.. code-block:: bash
   
    cat <fastq> | minerva_deconvolve -k 20 -w 40 -d 8 -a 50 --remove-stopwords > ebc_assignments.tsv
    
    minerva_deconvolve --help

Installation
------------

.. code-block:: bash
   
    git clone <url>   
    cd minerva_barcode_deconvolution
    python setup.py develop


Licence
-------

MIT License

Authors
-------

David Danko

`minerva_barcoded_read_deconvolution` was written by `David C. Danko	 <dcd3001@med.cornell.edu]>`_.
