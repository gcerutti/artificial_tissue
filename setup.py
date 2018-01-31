#!/usr/bin/env python
# -*- coding: utf-8 -*-

# {# pkglts, pysetup.kwds
# format setup arguments

from setuptools import setup, find_packages


short_descr = "A set of tools to generate virtual tissues with various possible computational representations"
readme = open('README.rst').read()
history = open('HISTORY.rst').read().replace('.. :changelog:', '')


# find version number in src/vplants/artificial_tissue/version.py
version = {}
with open("src/vplants/artificial_tissue/version.py") as fp:
    exec(fp.read(), version)


setup_kwds = dict(
    name='vplants.artificial_tissue',
    version=version["__version__"],
    description=short_descr,
    long_description=readme + '\n\n' + history,
    author="Guillaume Cerutti, Hadrien Oliveri, ",
    author_email="guillaume.cerutti@inria.fr, hadrien.oliveri@inria.fr, ",
    url='https://github.com/Guillaume Cerutti/artificial_tissue',
    license='cecill-c',
    zip_safe=False,

    packages=find_packages('src'),
    package_dir={'': 'src'},
    install_requires=[
        ],
    tests_require=[
        "mock",
        "nose",
        ],
    entry_points={},
    keywords='',
    test_suite='nose.collector',
)
# #}
# change setup_kwds below before the next pkglts tag

# do not change things below
# {# pkglts, pysetup.call
setup(**setup_kwds)
# #}
