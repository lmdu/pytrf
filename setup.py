import os
import glob
from setuptools import setup, Extension

extension = Extension('pytrf',
	sources = glob.glob('src/*.c'),
)

description = (
	"pytrf is a python package for finding tandem repeats from genomic sequences"
)

with open('README.rst') as fh:
	long_description = fh.read()

with open(os.path.join('src', 'version.h')) as fh:
	version = fh.read().split()[2].strip('"')

setup(
	name = 'pytrf',
	version = version,
	ext_modules = [extension],
	description = description,
	long_description = long_description,
	author = 'Lianming Du',
	author_email = 'adullb@qq.com',
	url = 'https://github.com/lmdu/pytrf',
	license = 'MIT',
	keywords = 'bioinformatics microsatellite tandem repeats',
	classifiers = [
			"Development Status :: 5 - Production/Stable",
			"Intended Audience :: Developers",
			"Intended Audience :: Education",
			"Intended Audience :: Science/Research",
			"Natural Language :: English",
			"License :: OSI Approved :: MIT License",
			"Programming Language :: C",
			"Programming Language :: Python :: 3.6",
			"Programming Language :: Python :: 3.7",
			"Programming Language :: Python :: 3.8",
			"Programming Language :: Python :: 3.9",
			"Programming Language :: Python :: 3.10",
			"Programming Language :: Python :: 3.11",
			"Programming Language :: Python :: 3.12",
			"Programming Language :: Python :: 3.13",
			"Operating System :: Microsoft :: Windows",
			"Operating System :: POSIX :: Linux",
			"Operating System :: Unix",
			"Topic :: Scientific/Engineering :: Bio-Informatics"
	],
	entry_points = {
		'console_scripts': ['pytrf = pytrfcli:main']
	},
	py_modules = ["pytrfcli"],
	install_requires = [
		"pyfastx>=2.1.0"
	],
	test_suite = "tests"
)