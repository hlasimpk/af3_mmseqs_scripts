from distutils.command.build import build
from distutils.util import convert_path
from setuptools import setup, find_packages
import sys
import numpy as np


class BuildCommand(build):
    user_options = build.user_options + [
        (
            "script-python-path=",
            None,
            "Path to Python interpreter to be included in the scripts",
        )
    ]

    def initialize_options(self):
        build.initialize_options(self)
        self.script_python_path = None

    def finalize_options(self):
        build.finalize_options(self)

    def run(self):
        global script_python_path
        script_python_path = self.script_python_path
        build.run(self)


def dependencies():
    with open("requirements.txt", "r") as f_in:
        deps = f_in.read().splitlines()
    return deps


def readme():
    with open("README.md", "r") as f_in:
        return f_in.read()


# ==============================================================
# Determine the Python executable
# ==============================================================
PYTHON_EXE = None
for arg in sys.argv:
    if arg[0:20] == "--script-python-path" and len(arg) == 20:
        option, value = arg, sys.argv[sys.argv.index(arg) + 1]
        PYTHON_EXE = value
    elif arg[0:20] == "--script-python-path" and arg[20] == "=":
        option, value = arg[:20], arg[21:]
        PYTHON_EXE = value

if not PYTHON_EXE:
    PYTHON_EXE = sys.executable

AUTHOR = "Adam Simpkin & Luc Elliott"
AUTHOR_EMAIL = "hlasimpk@liverpool.ac.uk"
DESCRIPTION = "Alphafold3 input processing tools"
DEPENDENCIES = dependencies()
LICENSE = "BSD License"
LONG_DESCRIPTION = readme()
PACKAGE_DIR = "af3_mmseqs2"
PACKAGE_NAME = "af3_mmseqs2"
PLATFORMS = ["Mac OS", "Windows", "Unix"]
URL = ""

PACKAGES = find_packages()


CLASSIFIERS = [
    "Development Status :: 4 - Beta",
    "Intended Audience :: Science/Research",
    "License :: OSI Approved :: BSD License",
    "Operating System :: OS Independent",
    "Programming Language :: Python",
    "Programming Language :: Python :: 3.6",
    "Programming Language :: Python :: 3.7",
    "Programming Language :: Python :: 3.8",
    "Programming Language :: Python :: 3.9",
    "Topic :: Scientific/Engineering :: Bio-Informatics",
]

setup(
    cmdclass={"build": BuildCommand},
    author=AUTHOR,
    author_email=AUTHOR_EMAIL,
    name=PACKAGE_NAME,
    description=DESCRIPTION,
    include_dirs=[np.get_include()],
    long_description=LONG_DESCRIPTION,
    license=LICENSE,
    url=URL,
    packages=PACKAGES,
    package_dir={PACKAGE_NAME: PACKAGE_DIR},
    data_files = ['data/config.ini'],
    platforms=PLATFORMS,
    classifiers=CLASSIFIERS,
    install_requires=DEPENDENCIES,
    include_package_data=True,
    zip_safe=False,
    entry_points={
        'console_scripts': [
            'alphafold3=af3_mmseqs2.alphafold3:main',
        ],
    },
)
