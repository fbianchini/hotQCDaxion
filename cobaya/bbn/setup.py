from setuptools import setup
import os

file_dir = os.path.abspath(os.path.dirname(__file__))
os.chdir(file_dir)

setup(name="bbn",
      version='1.0',
      description='BBN Cobaya likelihood package',
      zip_safe=False,  # set to false if you want to easily access bundled package data files
      packages=['bbn'],#, 'test_package.sub_module', 'test_package.tests'],
      package_data={'bbn': ['*.yaml', '*.bibtex', 'data/*',]},# 'data/**/*']},
      install_requires=['cobaya (>=2.0.5)'],
      # test_suite='test_package.tests',
      )

