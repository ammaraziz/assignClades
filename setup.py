from setuptools import setup

setup(name='CladeAssigner',
      version='0.1',
      description='Return HA clades for given input sequences',
      author='Ammar Aziz',
      author_email='ammar.aziz@mh.org.au',
      license='MIT',
      packages=['assignClades.py'],
      install_requires=['biopython', 'pandas'],
      zip_safe=False)