from setuptools import setup
from version import get_git_version

setup(name='dsa110-T2',
      version=get_git_version(),
      url='http://github.com/dsa110/dsa110-T2',
      packages=['T2'],
      zip_safe=False)
