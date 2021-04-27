import os
import glob
from setuptools import setup, Extension

extension = Extension('stria',
	sources = glob.glob('src/*.c'),
)

description = (
	"stria (short tandem repeat identification and analysis) is a "
	"python package for finding tandem repeats from DNA sequences"
)

with open('README.rst') as fh:
	long_description = fh.read()

with open(os.path.join('src', 'version.h')) as fh:
	version = fh.read().split()[2].strip('"')

setup(
	name = 'stria',
	version = version,
	ext_modules = [extension],
	description = description,
	long_description = long_description,
	author = 'Lianming Du',
	author_email = 'adullb@qq.com',
	url = 'https://github.com/lmdu/stria',
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
			"Programming Language :: Python :: 3.5",
			"Programming Language :: Python :: 3.6",
			"Programming Language :: Python :: 3.7",
			"Programming Language :: Python :: 3.8",
			"Programming Language :: Python :: 3.9",
			"Operating System :: Microsoft :: Windows",
			"Operating System :: POSIX :: Linux",
			"Operating System :: Unix",
			"Topic :: Scientific/Engineering :: Bio-Informatics"
	],
	entry_points = {
		'console_scripts': ['stria = striacli:main']
	},
	py_modules = ["striacli"],
	test_suite = "tests"
)