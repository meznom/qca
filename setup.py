#!/usr/bin/env python

# This is largely taken from the python-ecdsa project,
# https://github.com/warner/python-ecdsa.
import os, sys, subprocess, re
from distutils.core import setup, Command, Extension

VERSION_PY_FILENAME = 'qca/_version.py'
VERSION_PY = """\
# This file is originally generated from Git information by running 'setup.py
# version'. Distribution tarballs contain a pre-generated copy of this file.

__version__ = '{}'
"""

def update_version_py():
    """Update version file with version from Git.

    Constructs semver.org compatible version strings, provided that the Git tag
    is already a semver.org compatible version. The number of commits since the
    last tag and the abbreviated commit hash are appended as a pre-release
    version. An example: 2.0.0-25.gf52f551.
    """
    ver = 'unknown'
    try:
        p = subprocess.Popen(["git", "describe", "--dirty", "--always"],
                             stdout=subprocess.PIPE)
        stdout = p.communicate()[0]
        if p.returncode != 0:
            print('unable to run git')
        else:
            # Git describe returns one or more strings separated by a dash
            s = stdout.strip()
            fs = s.split('-')
            ver = fs[0]
            if len(fs) > 1:
                ver += '-' + '.'.join(fs[1:])
    except EnvironmentError:
        print('unable to run git')
    if os.path.exists(VERSION_PY_FILENAME) and ver == 'unknown':
        return
    f = open(VERSION_PY_FILENAME, "w")
    f.write(VERSION_PY.format(ver))
    f.close()
    print('set {} to "{}"'.format(VERSION_PY_FILENAME, ver))

def get_version():
    try:
        f = open(VERSION_PY_FILENAME)
    except EnvironmentError:
        return None
    for line in f.readlines():
        mo = re.match("__version__ = '([^']+)'", line)
        if mo:
            ver = mo.group(1)
            return ver
    return None

class Version(Command):
    description = 'update _version.py from Git repo'
    user_options = []
    boolean_options = []
    def initialize_options(self):
        pass
    def finalize_options(self):
        pass
    def run(self):
        update_version_py()
        print 'Version is now', get_version()

setup(name='qca',
      version=get_version(),
      description='QCA exact diagonalization.',
      author='Burkhard Ritter',
      author_email='burkhard@ualberta.ca',
      url='',
      packages=['qca','qca.test'],
      package_data={'qca': ['_qca.so']},
      cmdclass={ 'version': Version}
      )
