# Use setuptools, falling back on provided
try:
    from setuptools import setup, find_packages
except ImportError:
    import distribute_setup
    distribute_setup.use_setuptools()
    from setuptools import setup, find_packages

setup(name='anoisetools',
      version='0.2',
      author='Connor McCoy',
      author_email='cmccoy@fhcrc.org',
      packages=find_packages(),
      requires=['biopython (>=1.57)'],
      entry_points={
        'console_scripts': [
            'anoise = anoisetools.scripts.dispatch:main',
      ]},
      classifiers=[
        'License :: OSI Approved :: GNU General Public License (GPL)',
        'Development Status :: 3 - Alpha',
        'Programming Language :: Python :: 2.7',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
      ],
      license="GPL V3",
)
