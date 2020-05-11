from setuptools import setup
from os import path

VERSION = '0.6.2'

PACKAGES = [
    'py4pm',
]

# read the contents of your README file
this_directory = path.abspath(path.dirname(__file__))
with open(path.join(this_directory, 'README.md'), encoding='utf-8') as f:
        long_description = f.read()

setup(
    name='py4pm',
    version=VERSION,
    packages=PACKAGES,
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://gricad-gitlab.univ-grenoble-alpes.fr/webersa/py4pm',
    install_requires=['pandas', 'xlrd', 'matplotlib'],
    python_requires='>=3',
    author='Samuël Weber',
    author_email='samuel.weber@univ-grenoble-alpes.fr',
    project_urls={
        'Documentation': 'https://py4pm.readthedocs.io/en/latest/',
        'Source': 'https://gricad-gitlab.univ-grenoble-alpes.fr/webersa/py4pm',
	},
    classifiers=[
	# How mature is this project? Common values are
	#   3 - Alpha
	#   4 - Beta
	#   5 - Production/Stable
	'Development Status :: 3 - Alpha',

	# Indicate who your project is intended for
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Atmospheric Science',

	# Specify the Python versions you support here. In particular, ensure
	# that you indicate whether you support Python 2, Python 3 or both.
	'Programming Language :: Python :: 3',
	],
)
