
from setuptools import setup

# First the ParmedTools packages:
packages = ['mdoutanalyzer','pyrespgen']

modules = []

# Scripts
scripts = ['cpinutil.py', 'ceinutil.py', 'cpeinutil.py', 'charmmlipid2amber.py', 'mdout_analyzer.py',
           'finddgref.py', 'fitpkaeo.py', 'fixremdcouts.py', 'genremdinputs.py', 'softcore_setup.py',
           'py_resp.py', 'bar_pbsa.py','pyresp_gen.py']

if __name__ == '__main__':

    setup(name='AmberUtils',
          version='21.0',
          description='Various modules needed for AmberTools Python programs',
          author='Jason M. Swails, Ben Madej, Thomas T. Joseph, and Vinicius Wilian D. Cruzeiro',
          url='http://ambermd.org',
          license='GPL v2 or later',
          packages=packages,
          py_modules=modules,
          scripts=scripts)
