from distutils.core import setup

setup(name='gr50',
      version='0.1',
      description='Growth response metric tools',
      author='Jeremy Muhlich',
      author_email='jmuhlich@bitflood.org',
      url='https://github.com/sorgerlab/gr50_tools/',
      packages=['gr50'],
      scripts=['scripts/add_gr_column.py', 'scripts/compute_gr_metrics.py'],
      )
