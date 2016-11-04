#!/usr/bin/env python

from distutils.core import setup
import os

python_files = os.listdir('./pythonToolkit/')
moduleNames = []
for i in range(len(python_files)):
    if ".py" in python_files[i] and not "init" in python_files[i] and not ".swp" in python_files[i]:
	moduleNames.append('pythonToolkit.'+python_files[i].strip()[0:-3])

print moduleNames

setup(name='pythonToolkit',
      version='0.01',
      description='Assorted Python tools ',
      author='Christian Brinch',
      author_email='brinch@nbi.dk',
      py_modules=moduleNames)

