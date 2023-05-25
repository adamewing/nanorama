#!/usr/bin/env python

from setuptools import setup, find_packages

setup(
    name='ont-lrm',
    version='0.1.0',
    author='Adam Ewing',
    author_email='adam.ewing@gmail.com',
    description=("live read monitor"),
    license='MIT',
    url='https://github.com/adamewing/ont-lrm',
    download_url='https://github.com/adamewing/ont-lrm/archive/refs/tags/0.1.0.tar.gz',
    scripts=['lrm_mapper.py', 'lrm_viewer.py'],
    packages=find_packages(),
    install_requires = [
        'pysam',
        'pandas',
        'plotly',
        'dash',
        'bx-python',
        'numpy',
        'psutil',
    ],
    classifiers=[
        'Development Status :: 4 - Beta',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
    ],

)
