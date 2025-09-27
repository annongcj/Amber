# cd updateutils
# pytest -vs --cov=updateutils --cov-report=html .

import os
from os import path
import sys
from time import sleep 
import unittest
from mock import patch
sys.path.insert(0, os.path.abspath('..'))

from updateutils import main
from updateutils.main import _AMBER_VERSION, _AMBERTOOLS_VERSION
from updateutils.progressbar import AsciiProgressBar
from updateutils.upgrade import Upgrader
from updateutils.downloader import urlopen
from updateutils.patch import LocalPatch
from updateutils.patchlist import LocalPatchList


def test_main():
    main.MainUpdater()


def test_upgrade():
    _AMBERTOOLS_VERSION = 16
    _AMBER_VERSION = 16
    url = 'https://ambermd.org/upgrades/%d.%d/' % (_AMBERTOOLS_VERSION, _AMBER_VERSION)
    upgrade = Upgrader(url)
    upgrade.look_for_upgrade()
    upgrade.load_upgrade_info()


def test_progress_bar():
    progress = AsciiProgressBar(128)
    while not progress.finished: 
        sleep(0.02) 
        progress.update_add(8) 
    progress.reset() 
    while not progress.finished: 
        sleep(0.02) 
        progress.update_add(16)
    progress.reset() 
    while not progress.finished: 
        sleep(0.1) 
        progress.update_add(32)
