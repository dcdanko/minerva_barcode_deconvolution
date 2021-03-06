import setuptools

requirements = [
    'scikit-learn',
    'numpy',
    'pandas',
    'scipy',
    'click',
    'argparse',
]


setuptools.setup(
    name="minerva_deconvolve",
    version="1.1.3",
    url="https://github.com/dcdanko/minerva_barcode_deconvolution",
    author="David C. Danko",
    author_email="dcd3001@med.cornell.edu",
    description="Algorithm for deconvolving and clustering barcoded short reads",
    packages=setuptools.find_packages(),
    package_dir={'minerva': 'minerva'},
    entry_points={
        'console_scripts': [
            'minerva_deconvolve=minerva.deconvolution.deconvolve_barcodes:main',
            'minerva_enhance_kraken=minerva.kraken.enhance_kraken:main',
            'minerva_eval=minerva.eval.eval_deconvolution:main',
            'minerva_annotate=minerva.eval.annotate_fastq:main',
            'minerva_deconvolve_fastq=minerva.deconvolution.deconvolve_fasta:main',
        ]
    },
    install_requires=requirements,
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
    ],
)
