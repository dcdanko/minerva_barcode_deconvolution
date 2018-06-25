minerva_barcoded_read_deconvolution
===================================


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

Output
------

Minerva outputs a tsv file with three columns: read id, barcode, and cluster id. The deconvolved barcode for a read is a tuple of (barcode, cluster id).

Datasets
--------

The datasets used in the paper may be downloaded from AWS

`Dataset 1 <https://s3.us-east-2.amazonaws.com/minerva-datasets/10M.data1_atgctgaaq.fq.gz>`_


`Dataset 2 <https://s3.us-east-2.amazonaws.com/minerva-datasets/10M.data2_accctcct.fq.gz>`_

Licence
-------

MIT License

Authors
-------

David Danko

`minerva_barcoded_read_deconvolution` was written by `David C. Danko	 <dcd3001@med.cornell.edu]>`_.
