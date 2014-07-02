from setuptools import setup

setup(name='Seance',
      version='0.1',
      description='Reference-based phylogenetic analysis of 18S rRNA studies',
      author='Alan Medlar',
      author_email='alan.j.medlar@helsinki.fi',
      url='',
      license='GNU Public License ver3 ( https://www.gnu.org/licenses/gpl-3.0.html )',
      long_description='Evolutionary Analysis and Visualisation of Metagenomic Datasets',
      platforms=['*nix'],
      packages=['seance'],
      install_requires=[
          'biopython >= 1.6', 
          'dendropy >= 3.12',
          'cairocffi >= 0.5.4'
          ],
      scripts=['scripts/seance'],
     )

