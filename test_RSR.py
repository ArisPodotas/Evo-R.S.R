from cs50_final_project import Person, Population, compare, NCBI_parse, make_folders, cmd_line_input, make_figures, Sequence, Gene, Genome, Chromosome, Intergenic 
import pytest
import os
import random
import sys
from unittest.mock import patch

def main():
	test_sequence()
   
	test_gene()

	test_chromosome()

	test_genome()

	test_intragenic()

	test_Person_1()

	test_mut_1()

	test_compare_1()

	test_Population_1()

	test_make_folders()

	test_generations()

	test_cmd_line_input()

	test_cross()

	test_ncbi_parse()

	test_make_figures()

def test_sequence():
	sry = Sequence("ATATGCGC")
	assert sry.complement() == "TATACGCG" 
	assert sry.reverse() == "CGCGTATA"
	assert sry.reverse_complement() == "GCGCATAT"
	_ = Sequence("ATGC")
	assert _.point_mutations() != "ATGC"
	_ = Sequence("ATGC")
	assert _.insertions() != "ATGC" and len(_) == 5 
	_ = Sequence("ATGC")
	assert _.deletions() != "ATGC" and len(_) == 3
	_ = Sequence("ATGC")
	assert _.mutations() != "ATGC"
	_ = Sequence("A")
	assert _.transitions() == "G"
	_ = Sequence("ATGC")
	assert _.transversions() != "ATGC"
	_ = Sequence("ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC")
	yes = len(_)
	assert _.translocations() != "ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC" and len(_) == yes 
	_ = Sequence("ATGC")
	assert _.percent_contents() == [25, 25, 25, 25]
	assert _.adenine_percent() == 25
	assert _.thymine_percent() == 25 
	assert _.guanine_percent() == 25 
	assert _.cytocine_percent() == 25
	assert _.count_contents() == [1, 1, 1, 1]
	assert _.adenine_count() == 1
	assert _.guanine_count() == 1
	assert _.thymine_count() == 1
	assert _.cytocine_count() == 1
	p = Sequence("ATGTGGTAA")
	assert p.orf_finder(False, 0) == [["ATGTGGTAA", "MW"]]
	assert p.rearrange(1, 5, 3) == "AGTTGTGAA" 
	y = Sequence("AAAATG")

def test_gene():
	sry = Gene("ATGTAGTAGCCGAGCTGTTGTTGTGACACACACTGTGTGATGCACATATGCGAGTCTCTGCGAGCAGACGTATGGAGATAGGTGTGCGAGTGAGTAGAGCGCGCGTCGT", [[10, 15], [20, 25]], "test", "for the test")
	assert sry.name == "test"
	assert sry.description == "for the test"
	assert sry.seq == "ATGTAGTAGCCGAGCTGTTGTTGTGACACACACTGTGTGATGCACATATGCGAGTCTCTGCGAGCAGACGTATGGAGATAGGTGTGCGAGTGAGTAGAGCGCGCGTCGT"
	assert sry.exons == [[10, 15], [20, 25]]
	assert sry.introns == [15, 20]
	assert sry.mature() == ["CGAGC", "TTGTG"]
	assert sry.intron_seq() == ["TGTTG"]

def test_intragenic():
	ran = Intergenic("ATGTGTCC", "YES")
	assert ran.name == "YES"
	assert ran.length == 8
	assert ran.seq == "ATGTGTCC" 

def test_chromosome():
	chrom = Chromosome("ATGCGTGTGGTGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG", [], 17, "")

def test_genome():
	gen = Genome("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAATTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC", [])

def test_Person_1():
	p0 = Person(sequence = random.choices(population = ["A", "T", "G", "C"], k = len(NCBI_parse("ATTTGCGTAG"))))
	with pytest.raises(TypeError):
		p0 = Person(sequence = 2)
		po = Person(ID = "7")
		po = Person(ID = [1, 2, 2])
		po = Person(ID = -2)
		po = Person(ID = complex(3))
	with pytest.raises(ValueError):
		p0 = Person(sequence = "boy")
		p0 = Person(genome = "2345432")
		p0 = Person(genome = 33)
	p = Person(sequence = random.choices(population = ["A", "T", "G", "C"], k = len(NCBI_parse("ATTTGCGTAG"))))
	assert p.genome == 10
	per = Person(sequence = random.choices(population = ["A", "T", "G", "C"], k = len(NCBI_parse("ATTTGCGTAG"))))
	assert per.seq != None
	assert per.seq != ""

def test_mut_1():
	pm = Person(sequence = random.choices(population = ["A", "T", "G", "C"], k = len(NCBI_parse("ATTTGCGTAG"))))
	assert pm.seq != pm.mutate()
	with pytest.raises(ValueError):
		pm.mutate(rate = -3)

# def test_translocate_1():
#	 p = Person(reference_sequence = NCBI_parse("ATTTGCGTAG"), sequence = random.choices(population = ["A", "T", "G", "C"], k = len(NCBI_parse("ATTTGCGTAG"))))
#	 assert p.allele1 == p.allele2

# def test_translocate_2():
#	 p = Person(reference_sequence = NCBI_parse("ATTTGCGTAG"), sequence = random.choices(population = ["A", "T", "G", "C"], k = len(NCBI_parse("ATTTGCGTAG"))))
#	 sum = p.allele1 + p.allele2
#	 p.translocate()
#	 assert sum != p.allele2 + p.allele1

def test_compare_1():
	Person_compare = Person(sequence = random.choices(population = ["A", "T", "G", "C"], k = len(NCBI_parse("ATTTGCGTAG"))))
	assert compare(Person_compare.seq, Person_compare.seq) == 100
	before = Person_compare.seq
	Person_compare.mutate()
	after = Person_compare.seq
	assert 0 <= compare(before, after) <= 100

def test_Population_1():
	athens = Population()
	assert athens.log != []
	malta = Population()
	assert str(malta) != ""
	assert malta.reference.seq == "AAAGGTACGCGCGCCGGCGCGTATAGCTTTAGTCGTGGACGCTAGCTAGCTGGTAGCGACAGGCGAGAAATGCTAGCATCGAGCATGCAGCGTTC"
	with pytest.raises(ValueError):
		czechia = Population(generations = -3)
		andorra = Population(generations = "hgj")

# def test_Population_trees_1():
#	 athens = Population()
#	 assert " " in athens.tree(300, 300)

def test_make_folders():
	output_path = "."
	new_folder = "new_folder"
	make_folders(output_path, new_folder)
	assert os.path.exists(new_folder) == True

def test_generations():
	athens = Population(size = 100, generations = 10)
	result = athens.generations()
	assert isinstance(result[0], list) == True
	assert isinstance(result[1], list) == True
	assert isinstance(result[2], str) == True
	for i in range(len(athens.slog)):
		assert athens.slog[i] >= 0
	for n in range(len(athens.average)):
		assert isinstance(athens.average[n], float) == True
		assert 0.00 <= athens.average[n] <= 100.00

def test_cmd_line_input():
	assert cmd_line_input()['-o'] == "."
	assert cmd_line_input()['-acc'] == 'AAAGGTACGCGCGCCGGCGCGTATAGCTTTAGTCGTGGACGCTAGCTAGCTGGTAGCGACAGGCGAGAAATGCTAGCATCGAGCATGCAGCGTTC'
	assert cmd_line_input()['-gen'] == 100
	assert cmd_line_input()['-size'] == 100
	assert cmd_line_input()['-rate'] == 0.5
	assert cmd_line_input()["-v"] == False
	assert cmd_line_input()["-q"] == True
	# assert cmd_line_input()["-speed"] == 100
	# assert cmd_line_input()["-regression"] == False
	assert cmd_line_input()["-log"] == True
	assert cmd_line_input()["-drift"] == 10 
	assert cmd_line_input()["-criteria"] == 10 
	with patch.object(sys, 'argv', ['cs50_final_project.py', '-o', '.', '-acc', 'AAAT', '-gen', '5', '-size', '100', '-rate', '0.5', "-v", "True", "-log", "True", "-drift", "10", "-criteria", "10"]):
		assert cmd_line_input()['-acc'] == 'AAAT'
		assert cmd_line_input()['-gen'] == 5
		assert cmd_line_input()['-size'] == 100
		assert cmd_line_input()['-rate'] == 0.5
		assert cmd_line_input()["-v"] == True
		assert cmd_line_input()["-q"] == False
		# assert cmd_line_input()["-speed"] == 99
		# assert cmd_line_input()["-regression"] == True 
		assert cmd_line_input()["-log"] == True
		assert cmd_line_input()["-drift"] == 10 
		assert cmd_line_input()["-criteria"] == 10 
	with pytest.raises(SystemExit):
		with patch.object(sys, 'argv', ['cs50_final_project.py', "-h"]):
			cmd_line_input()
	with pytest.raises(SystemExit):
		with patch.object(sys, 'argv', ['cs50_final_project.py', '-acc', 'ATGGGG', "-w"]):
			cmd_line_input()

def test_cross():
	athens = Population()
	athens.cross(athens.sequences)
	assert athens.log != None

def test_ncbi_parse():
	assert NCBI_parse('LM378763.1') == "TAGCTTATCAGACTGATGTTGA"
	with pytest.raises(ValueError):
		NCBI_parse('hmnsou.1')
		NCBI_parse("ATGCTAGCTTGATCTGTPPPPP")
		NCBI_parse("ATTATATATATA88888888888")
		NCBI_parse("ATGCTAGCTTGATCTGT.1")

def test_make_figures():
	athens = Population(generations = 10)
	generations_output = athens.generations()
	assert make_figures(generations_output[0], generations_output[1]) == True

def test_Person():
	aris = Person()
	assert isinstance(aris, object) == True
	assert aris.seq != aris.mutate()	

if __name__ == "__main__":
	main()
