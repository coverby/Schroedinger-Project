from distutils.core import setup

with open('README.md') as f:
      long_description = ''.join(f.readlines())

setup(name='Langevin477CO',
      version='.1',
      description='CHE 477 Langevin Simulator',
      long_description=long_description,
      author='Clyde Overby',
      packages=['Langevin477CO'],
      entry_points=
      {
            'console_scripts': ["langevin=Langevin477CO.lang:start"]
      }
     )