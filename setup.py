from setuptools import setup, find_packages

setup(
    name='gibby_motifFinding',  # Replace with your package name
    version='0.9',  # Update to a new version number
    packages=find_packages(),
    entry_points={
        'console_scripts': [
            'gibby_run = mypackage.gibby_run:main'
        ]
    },
    install_requires=[
        'biopython',
        'numpy',
        'seqlogo'
    ],
    include_package_data=True,
    description='Process genomic data to extract peak sequences and compute PWMs.',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url='https://github.com/kairitanaka/CSE_185_finalProject',  # Replace with your GitHub URL
    author='Joe Hwang, Kairi Tanaka',  # Replace with your name
    author_email='j8hwang@ucsd.edu, ktanaka@ucsd.edu',  # Replace with your email
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.6',
)


