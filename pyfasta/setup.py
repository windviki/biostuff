from setuptools import setup, find_packages


version = '0.3.1'

setup(name='pyfasta',
      version=version,
      description="fast, memory-efficient, pythonic access to fasta sequence files",
      url="http://bitbucket.org/brentp/biostuff/",
      long_description=open('README.txt').read(),
      classifiers=["Topic :: Scientific/Engineering :: Bio-Informatics"],
      keywords='bioinformatics blast fasta',
      author='brentp',
      author_email='bpederse@gmail.com',
      license='MIT',
      packages=find_packages(exclude=['ez_setup', 'examples', 'tests']),
      include_package_data=True,
      test_suite='nose.collector',
      zip_safe=False,
      install_requires=[
          # -*- Extra requirements: -*-
      ],
      scripts=[],
      entry_points={
      'console_scripts': ['pyfasta = pyfasta:main']
      },
      )
