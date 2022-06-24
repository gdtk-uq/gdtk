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
    print("Creating directory: ", dir)
    copytree('blank', dir)
    chdir(dir)
    run('sed -i s/#OVERRIDE_CORE_COUNT/np={}/g makefile'.format(np).split(), check=True)
    chdir('../')

