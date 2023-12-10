# setup.py
from setuptools import setup, find_packages

setup(
    name='Astra',
    version='0.1',
    packages=find_packages(),
    install_requires=[
        'pyhmmer',
        'pandas',
        'tqdm',
        'requests',
	    'asyncio'
    ],
        entry_points={
            'console_scripts': ['Astra=astra.main:main'],
        },
)
