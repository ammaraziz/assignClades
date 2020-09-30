from setuptools import setup

setup(name='CladeAssigner',
      version='0.1',
      description='Return clades for given input sequences',
      author='Ammar Aziz',
      author_email='jrj.healey@gmail.com',
      license='MIT',
      packages=['prymer'],
      install_requires=['biopython', 'seqanpy', 'pandas'],
      zip_safe=False)