import setuptools

requirements = [
    'Click>=6.0',
    'scikit-learn==0.19.1',
    'numpy==1.13.3',
    'pandas==0.19.2'
    # TODO: put package requirements here
]


setuptools.setup(
    name="minerva_barcoded_read_deconvolution",
    version="0.1.0",
    url="https://github.com/dcdanko/minerva_barcode_deconvolution",

    author="David C. Danko",
    author_email="dcd3001@med.cornell.edu",

    description="Toolsets for deconvolving and clustering barcoded short reads",
    long_description=open('README.rst').read(),

    packages=setuptools.find_packages(),
    package_dir={'minerva': 'minerva'},
    
    ext_modules = [setuptools.Extension("cseqs", ["cext_minerva/_cseqs.c","cext_minerva/cseqs.c"])],
    
    entry_points = {
        'console_scripts': [
            'minerva_deconvolve=minerva.deconvolution.deconvolve_barcodes:main',
            'minerva_enhance_kraken=minerva.kraken.enhance_kraken:main',
            'minerva_eval=minerva.eval.eval_deconvolution:main',
            'minerva_annotate=minerva.eval.annotate_fastq:main',
        ]
    },
    
    install_requires=requirements,

    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
    ],
)
