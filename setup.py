from distutils.core import setup

with open('README.md') as f:
      long_description = ''.join(f.readlines())

setup(name='Schro1D',
      version='.01',
      description='CHE 477 1D Numerical Schroedinger Solver',
      long_description=long_description,
      author='Clyde Overby',
      packages=['Schro1D'],
      entry_points=
      {
            'console_scripts': ["schro=Schro1D.lang:start"]
      }
     )