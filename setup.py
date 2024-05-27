from setuptools import setup, find_packages

setup(
    name='gibby',
    version='0.1.0',
    packages=find_packages(),
    install_requires=[
        'biopython',
        'numpy',
        'seqlogo',
        'tqdm'
    ],
    entry_points={
        'console_scripts': [
            'gibby = gibby.gibby:main'
        ]
    },
    author='Kairi Tanaka, Joseph Hwang',
    author_email='ktanaka@ucsd.edu',
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



