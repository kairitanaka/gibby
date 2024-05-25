from setuptools import setup, find_packages

setup(
    name='your_package_name',
    version='0.1.0',
    packages=find_packages(),
    install_requires=[
        'biopython',
        'numpy',
        'seqlogo'
    ],
    entry_points={
        'console_scripts': [
            'gibby_motifFinding = your_package.gibby_motifFinding:main'
        ]
    },
    author='Your Name',
    author_email='your.email@example.com',
    description='A package for genomic motif finding.',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url='https://github.com/kairitanaka/CSE_185_finalProject',
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.6',
)
