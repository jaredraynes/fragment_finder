import pytest
import frag_finder

#to run pytest from the command line be in the directory with your program and then:> python -m pytest

vclh_seq = 'MNHKVHMHHHHHHADEQEEKAKVRTELIQELAQGLGGIEKKNFPTLGDEDLDHTYMTKLLTYLQEREQAENSWRKRLLKGIQDHALDLVPRGSPGLPGPRGEQGPTGPTGPAGPRGLQGLQGLQGERGEQGPTGPAGPRGLQGERGEQGPTGLAGKAGEAGAKGETGPAGPQGPRGEQGPQGLPGKDGEAGAQGRPGKRGKQGQKGEKGEPGTQGAKGDRGETGPVGPRGERGEAGPAGKDGERGPVGPAGKDGQNGQDGLPGKDGKDGQNGKDGLPGKDGKDGQNGKDGLPGKDGKDGQDGKDGLPGKDGKDGLPGKDGKDGQPGKPGKY'

def test_mass_cal():
	mass = frag_finder.mass_cal(vclh_seq)
	expected = 33766.85
	assert mass == expected

def test_fragments():
	fragments = frag_finder.fragments(vclh_seq)
	expected = ""
	assert fragments == expected

