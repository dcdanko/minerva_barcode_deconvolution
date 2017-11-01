import setuptools

requirements = [
    'Click>=6.0',
    # TODO: put package requirements here
]


setuptools.setup(
    name="minerva_barcoded_read_deconvolution",
    version="0.1.0",
    url="https://github.com/borntyping/cookiecutter-pypackage-minimal",

    author="David C. Danko	",
    author_email="dcd3001@med.cornell.edu",

    description="Toolsets for deconvolving and clustering barcoded short reads",
    long_description=open('README.rst').read(),

    packages=['minerva'],
    package_dir={'minerva': 'minerva'},
    
    ext_modules = [setuptools.Extension("cseqs", ["cext_minerva/_cseqs.c","cext_minerva/cseqs.c"])],
    
    entry_points = {
        'console_scripts': [
            'minerva_deconvolve=minerva.deconvolution.deconvolve_barcodes:main',
            'minerva_enhance_kraken=minerva.kraken.enhance_kraken:main'
            ]
    },
    
    install_requires=requirements,

    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Programming Language :: Python',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
    ],
)