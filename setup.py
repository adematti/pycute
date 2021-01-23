import os
import sys
from setuptools import setup
from setuptools.command.install import install
from distutils.command.build import build as _build
from setuptools.command.bdist_egg import bdist_egg as _bdist_egg
from subprocess import call

setup_keywords = dict()
setup_keywords['name'] = 'pycute'
setup_keywords['description'] = 'Correlation utility'
setup_keywords['author'] = 'Arnaud de Mattia'
setup_keywords['author_email'] = ''
setup_keywords['license'] = 'GPL3'
setup_keywords['url'] = 'https://github.com/adematti/pycute'
sys.path.insert(0,os.path.abspath('pycute/'))
import version
setup_keywords['version'] = version.__version__
setup_keywords['install_requires'] = ['numpy']
setup_keywords['packages'] = ['pycute']
setup_keywords['package_dir'] = {'pycute':'pycute'}

class build(_build):

    def run(self):
        super(build,self).run()
        def compile():
            call('make')
        self.execute(compile,[],'Compiling pycute')
        self.move_file('pycute/cute.so',os.path.join(self.build_lib,'pycute'))

class bdist_egg(_bdist_egg):
    def run(self):
        self.run_command('build')
        _bdist_egg.run(self)

setup_keywords['cmdclass'] = {'build':build,'bdist_egg':bdist_egg}

setup(**setup_keywords)
