import glob
from setuptools import setup, Extension

extension = Extension('stripy',
	sources = glob.glob('src/*.c'),
)

description = (
	"strit (short tandem repeat identification tool) is a python "
	"module for finding tandem repeats from DNA sequences"
)

with open('README.rst') as fh:
	long_description = fh.read()

version = '0.0.1'

setup(
	name = 'stripy',
	version = version,
	ext_modules = [extension],
	description = description,
	long_description = long_description,
	author = 'Lianming Du',
	author_email = 'adullb@qq.com',
	url = 'https://github.com/lmdu/stripy',
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
		'console_scripts': ['stripy = stripycli:main']
	},
	py_modules = ["stripycli"],
	test_suite = "tests"
)