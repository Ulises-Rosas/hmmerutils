import setuptools
from distutils.core import setup

setup(
    name = 'hmmerutils',
    version = '0.1',
    author = 'Ulises Rosas',
    packages = ['hmmerutils'],
    entry_points={
        'console_scripts': [
            'hmmerparser  = hmmerutils.hmmerutils_core:main'
            ]
    },
    classifiers = [
        'Programming Language :: Python :: 3'
        ]
)

