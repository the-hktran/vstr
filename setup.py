#!/usr/bin/env python

"""The setup script."""

from setuptools import setup, find_packages

with open('README.md') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()

requirements = [ ]

test_requirements = [ ]

setup(
    author="Henry Khoa Tran",
    author_email='hkt2112@columbia.edu',
    python_requires='>=3.6',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
    ],
    description="Vibrational structure theory methods based on heat bath CI",
    entry_points={
        'console_scripts': [
            'vstr=vstr.cli:main',
        ],
    },
    install_requires=requirements,
    license="MIT license",
    long_description=readme + '\n\n' + history,
    include_package_data=True,
    keywords='vstr',
    name='vstr',
    packages=find_packages(include=['vstr', 'vstr.*']),
    test_suite='tests',
    tests_require=test_requirements,
    url='https://github.com/the-hktran/VibrationalStructure',
    version='0.1.0',
    zip_safe=False,
)
