import argparse
import numpy
import matplotlib
import lxml
import pandas
import pyteomics


#A program to identify peptide fragments and their location within a protein 
#using mass spectrometry data

# https://pyteomics.readthedocs.io/en/latest/index.html

# potential code to use:

>>> from pyteomics import parser
>>> parser.cleave('AKAKBK', parser.expasy_rules['trypsin'], 0)
{'AK', 'BK'}

from pyteomics import parser
forms = parser.isoforms('PEPTIDE', variable_mods={'p': ['T'], 'ox': ['P']})
for seq in forms: print seq

>>> mass.calculate_mass(sequence='PEPTIDE')
799.359964027207

>>> from pyteomics import parser
>>> ps = parser.parse('PEPTIDE', show_unmodified_termini=True)
>>> mass.calculate_mass(parsed_sequence=ps)
799.359964027207

>>> from pyteomics import mass
>>> mass.calculate_mass(sequence='PEPTIDE', ion_type='M', charge=2)
400.6872584803735

>>> mass.calculate_mass(sequence='PEP', ion_type='b', charge=1)
324.15539725264904

>>> mass.calculate_mass(sequence='TIDE', ion_type='y', charge=1)
477.219119708098




