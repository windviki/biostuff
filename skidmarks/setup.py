from distutils.core import setup
from setuptools import setup

version="0.0.2"

setup(name='skidmarks',
      version=version,
     description="find runs (non-randomness) in sequences",
      url="http://code.google.com/p/bpbio/",
      long_description=open('README.txt').read(),
      classifiers=["Topic :: Scientific/Engineering :: Bio-Informatics"],
      keywords='bioinformatics sequence randomness test',
      author='brentp',
      author_email='bpederse@gmail.com',
      license='MIT',
      packages=['.'],
      test_suite='skidmarks.test_suite',
      zip_safe=False,
      install_requires=[
          # -*- Extra requirements: -*-
      ],
      scripts=[],
      entry_points={
      },
)


