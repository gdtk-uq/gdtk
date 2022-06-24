"""
Automatic run script for Eilmer parallel scaling tests.
"""

from shutil import copytree
from subprocess import run
from os import chdir
import configparser

config = configparser.ConfigParser()
config.read('parameters.txt')
nps = config['parameters']['number_of_processors'].split(',')

for np in nps:
    dir = 'test.{}'.format(str(np).zfill(2))
    chdir(dir)
    run('make partition'.split(), check=True)
    run('make prep'.split(), check=True)
    run('make run'.split(), check=True)
    chdir('../')

