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
    version="1.0.0",
    url="https://github.com/dcdanko/minerva_barcode_deconvolution",
    author="David C. Danko",
    author_email="dcd3001@med.cornell.edu",
    description="Algorithm for deconvolving and clustering barcoded short reads",
    packages=setuptools.find_packages(),
    package_dir={'minerva': 'minerva'},
    ext_modules=[
        setuptools.Extension("cseqs", ["cext_minerva/_cseqs.c", "cext_minerva/cseqs.c"])
    ],
    entry_points={
        'console_scripts': [
            'minerva_deconvolve=minerva.deconvolution.deconvolve_barcodes:main',
            'minerva_enhance_kraken=minerva.kraken.enhance_kraken:main',
            'minerva_eval=minerva.eval.eval_deconvolution:main',
            'minerva_annotate=minerva.eval.annotate_fastq:main',
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
