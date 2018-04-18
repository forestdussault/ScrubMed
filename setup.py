#!/usr/bin/env python

from setuptools import setup, find_packages

setup(
    name="ScrubMed",
    version="0.1.3",
    packages=find_packages(),
    author="Forest Dussault",
    author_email="forest.dussault@inspection.gc.ca",
    url="https://github.com/forestdussault/ScrubMed",
    scripts=['scrubmed.py'],
    install_requires=['click', 'pandas', 'bs4', 'BioPython', 'lxml']
)
