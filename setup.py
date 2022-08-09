from setuptools import setup
#from version import get_git_version

setup(name='GReX-T2',
      version='0.1.0',
      url='https://github.com/GReX-Telescope/GReX-T2.git',
      requirements=['seaborn', 'astropy', 'hdbscan', 'progress'],
      packages=['T2'],
      zip_safe=False)
