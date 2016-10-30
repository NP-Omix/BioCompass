#!/usr/bin/env python
# -*- coding: utf-8 -*-

from setuptools import setup

with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()

with open('requirements.txt') as requirements_file:
        requirements = requirements_file.read()

test_requirements = [
    # TODO: put package test requirements here
]

setup(
    name='BioCompass',
    version='0.6.0',
    description="Python package for gene clustering",
    long_description=readme + '\n\n' + history,
    author="Tiago Leão, Gui Castelão",
    author_email='tferreir@ucsd.edu, guilherme@castelao.net',
    url='https://github.com/NP-Omix/BioCompass',
    packages=[
        'BioCompass',
    ],
    package_dir={'BioCompass':
                 'BioCompass'},
    entry_points={
        'console_scripts': [
            'BioCompass=BioCompass.cli:main'
        ]
    },
    include_package_data=True,
    install_requires=requirements,
    license="BSD license",
    zip_safe=False,
    keywords='BioCompass',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: BSD License',
        'Natural Language :: English',
        "Programming Language :: Python :: 2",
        'Programming Language :: Python :: 2.6',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
    ],
)
