# Project name: R.S.R.
# Collaborators Rafail Adam, Aris Podotas
# Country: Greece
# City: Athens
# Date: 30/8/2023

import random
import re
import os
import math
from Bio import Entrez
from Bio import SeqIO
import sys
import matplotlib.pyplot as plt
import statistics as stat
import time
from tkinter import *
from tkinter import ttk
import bigtree as bt
import numpy as np

class Sequence:
	"""This class is for representing sequences of DNA or RNA."""
	def __init__(self, sequence: str|list):
		if isinstance(sequence, str):
			self.seq = sequence
			self.seqls = list(self.seq)
		elif isinstance(sequence, list):
			self.seq = ''.join(sequence)
			self.seqls = sequence
		else:
			raise TypeError(f"The sequence must be a string or list.\nExpecter type (str, list).\nGot {type(sequence)}.\n")
		if re.search(r"[AGCT]+", self.seq):
			self.dna = True
			self.rna = False
		elif re.search(r"[AUGC]+", self.seq):
			self.rna = True
			self.dna = False
		else:
			raise ValueError(f"False sequence character in sequence.\n")
		
	def __str__(self):
		"""Returns the sequence of the instance."""
		return self.seq

	def __len__(self):
		"""Returns the sequence lenght."""
		return len(self.seq)
	
	def __add__(self, other):
		if self.dna and other.dna:
			seq = self.seq + other.seq
		elif self.rna and other.rna:
			seq = self.seq + other.seq
		else:
			raise ReferenceError("Cannot concatinate instances of rna and dna.\n")
		return seq
	
	def complement(self):
		"""Returns the complementary sequence.
		Does not alter the original sequence."""
		self.comp_seq = ""
		if self.dna == True:
			for nuc in self.seq:
				if nuc == "A":
					self.comp_seq += "T"
				if nuc == "T":
					self.comp_seq += "A"
				if nuc == "G":
					self.comp_seq += "C"
				if nuc == "C":
					self.comp_seq += "G"
		elif self.rna == True:
			for nuc in self.seq:
				if nuc == "A":
					self.comp_seq += "U"
				if nuc == "U":
					self.comp_seq += "A"
				if nuc == "G":
					self.comp_seq += "C"
				if nuc == "C":
					self.comp_seq += "G"
		return self.comp_seq
	
	def reverse(self):
		"""Returns the reversed sequence.
		Does not alter the original sequence."""
		self.rev_seq = self.seq[::-1]
		return self.rev_seq

	def reverse_complement(self):
		"""Returns the reversed and complementary sequence.
		Does not alter the original sequence"""
		self.rc_seq = Sequence(sequence = self.complement()).reverse()
		return self.rc_seq

	def insertions(self, ammount: int = 1):
		"""Returns a mutated version of the sequence where an insertion has happened.
		Alters original sequence"""
		if not isinstance(ammount, int):
			raise TypeError(f"Please specify a positive ammount of mutations.\n")
		if ammount <= 0:
			raise ValueError(f"The ammount of mutations must be positive and whole.\n")
		for _ in range(ammount):
			seed = random.randint(0, len(self.seq) - 1)
			self.seqls[seed:seed + 1:1] = f'{self.seqls[seed]}{random.choice(["A", "T", "G", "C"])}'
			self.seq = ''.join(self.seqls)
		return self.seq

	def deletions(self, ammount: int = 1):
		"""Returns a mutated version of the sequence where a deletion has happened.
		Alters the original sequence"""
		if not isinstance(ammount, int):
			raise TypeError(f"Please specify a positive ammount of mutations.\n")
		if ammount <= 0:
			raise ValueError(f"The ammount of mutations must be positive and whole.\n")
		for _ in range(ammount):
			seed = random.randint(0, len(self.seq) - 1)
			self.seqls[seed] = ""
			self.seq = ''.join(self.seqls)
		return self.seq

	def mutations(self, ammount: int = 1, mutation_weights: dict = {"point mutations": 20, "transversions": 40, "transitions": 60, "insertions": 80, "deletions": 100, "translocations": 120}):
		"""Returns a sequence with random types of mutations at a specific rate.
		Alters the original sequence."""
		if not isinstance(ammount, int):
			raise TypeError(f"Please specify a positive ammount of mutations.\n")
		if ammount <= 0:
			raise ValueError(f"The ammount of mutations must be positive and whole.\n")
		if not isinstance(mutation_weights, dict):
			raise TypeError(f"Please use a dictionary to specify mutation weights.\n")
		for key in mutation_weights:
			if key not in ["point mutations", "transversions", "transitions", "insertions", "deletions", "translocations"]:
				raise ValueError(f"Wrong type of mutation in mutation weights: {key}.\n")
		for value in mutation_weights.values():
			if value < 0:
				raise ValueError(f"Invalid weight for mutation: {value}.\n")
		modified_mutation_weights = [mutation_weights["point mutations"]]
		modified_mutation_weights.append(mutation_weights["transversions"])
		modified_mutation_weights.append(mutation_weights["transitions"])
		modified_mutation_weights.append(mutation_weights["insertions"])
		modified_mutation_weights.append(mutation_weights["deletions"])
		modified_mutation_weights.append(mutation_weights["translocations"])
		for _ in range(ammount):
			mut = random.choices([self.point_mutations,
						 self.transversions,
						 self.transitions,
						 self.insertions,
						 self.deletions,
						 self.translocations], cum_weights = modified_mutation_weights)
			mut[0](ammount = 1)
		return self.seq

	def point_mutations(self, ammount: int = 1):
		"""Returns the sequence with one point mutation.
		Alters the original sequence."""
		if not isinstance(ammount, int):
			raise TypeError(f"Please specify a positive ammount of mutations.\n")
		if ammount <= 0:
			raise ValueError(f"The ammount of mutations must be positive and whole.\n")
		ls = list(self.seq)
		for _ in range(ammount):
			mut_position = random.randint(0, len(ls) - 1)
			nuc_in_pos = self.seq[mut_position]
			if nuc_in_pos == "A":
				mut_options = ["T", "G", "C"]
			elif nuc_in_pos == "T":
				mut_options = ["A", "G", "C"]
			elif nuc_in_pos == "G":
				mut_options = ["A", "T", "C"]
			elif nuc_in_pos == "C":
				mut_options = ["A", "T", "G"]
			ls[mut_position] = random.choice(mut_options)
			self.seq = ''.join(ls)
		return self.seq
	
	def transversions(self, ammount: int = 1):
		"""Returns the sequence having done a transversion mutation.
		Alters the original sequence.
		A transversion is a mutation that changes a pyrimidine to a purine and reverse.
		Essentially turns A <-> T, C. G <-> T, C."""
		if not isinstance(ammount, int):
			raise TypeError(f"Please specify a positive ammount of mutations.\n")
		if ammount <= 0:
			raise ValueError(f"The ammount of mutations must be positive and whole.\n")
		ls = list(self.seq)
		for _ in range(ammount):
			mut_position = random.randint(0, len(self.seqls) - 1)
			nuc_in_pos = self.seq[mut_position]
			if nuc_in_pos == "A":
				mut_options = ["T", "C"]
			elif nuc_in_pos == "T":
				mut_options = ["G", "A"]
			elif nuc_in_pos == "G":
				mut_options = ["T", "C"]
			elif nuc_in_pos == "C":
				mut_options = ["G", "A"] 
			ls[mut_position] = random.choice(mut_options)
			self.seq = ''.join(ls)
		return self.seq
	
	def transitions(self, ammount: int = 1):
		"""Returns a sequence with a transition mutatuin.
		Alters the original sequence.
		A transition mutation is a mutation that changes a pyrimidine to a different pyrimidine or a purine to a different purine.
		Essentially A <-> G, C <-> T.""" 
		if not isinstance(ammount, int):
			raise TypeError(f"Please specify a positive ammount of mutations.\n")
		if ammount <= 0:
			raise ValueError(f"The ammount of mutations must be positive and whole.\n")
		ls = list(self.seq)
		for _ in range(ammount):
			mut_position = random.randint(0, len(self.seqls) - 1)
			nuc_in_pos = self.seq[mut_position]
			if nuc_in_pos == "A":
				mut_options = ["G"]
			elif nuc_in_pos == "T":
				mut_options = ["C"]
			elif nuc_in_pos == "G":
				mut_options = ["A"]
			elif nuc_in_pos == "C":
				mut_options =  ["T"]
			ls[mut_position] = random.choice(mut_options)
			self.seq = ''.join(ls)
		return self.seq
	
	def rearrange(self, start: int|str, stop: int|str, place: int):
		"""Returns a rearranged sequence."""
		if not isinstance(start, int|str):
			raise TypeError(f"Start should be a place in the sequnce. Given {start}.\n")
		if not isinstance(stop, int|str):
			raise TypeError(f"Stop should be a place in the sequence. Given {stop}.\n")
		if start >= stop or start == stop:
			raise ValueError(f"Invalid values of start and stop. With given (start > stop) the segment will be null. Given start: {start}, stop: {stop}")
		segment = self.seq[start:stop:1]
		self.seq = self.seq[:start:1] + self.seq[stop::1]
		self.seq = self.seq[:place:1] + segment + self.seq[place::1]
		return self.seq

	def translocations(self, ammount: int = 1):
		"""Returns a mutated version of the sequence where there has been a translocation of a segment.
		Alters the original sequence."""
		if not isinstance(ammount, int):
			raise TypeError(f"Please specify a positive ammount of mutations.\n")
		if ammount <= 0:
			raise ValueError(f"The ammount of mutations must be positive and whole.\n")
		for _ in range(ammount):
			seed_start = random.randint(0, len(self.seq) - 2)
			seed_end = random.randint(1, len(self.seq) - 1)
			if seed_start > seed_end:
				seed_start = seed_start - seed_end
				seed_end += seed_start
				seed_start = seed_end - seed_start 
			elif seed_start == seed_end:
				seed_end = random.randint(seed_start + 1, len(self.seq) - 1)
			segment = self.seq[seed_start:seed_end:1]
			self.seq = self.seq[:seed_start:1] + self.seq[seed_end::1]
			position = random.randint(0, len(self.seq))
			self.seq = self.seq[:position:1] + segment + self.seq[position::1]
		return self.seq

	def orf_finder(self, nested: bool = True, min_size: int = 0, start_codons: list = ["ATG"], stop_codons: list = ["TAA", "TAG", "TGA"], codon_size: int = 3, genetic_code: dict = {"TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L", "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S", "TAT": "Y", "TAC": "Y", "TAA": "", "TAG": "", "TGT": "C", "TGC": "C", "TGA": "", "TGG": "W", "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L", "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P", "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q", "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R", "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M", "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T", "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K", "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R","GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V", "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A", "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E", "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G"}):
		"""Returns all open reading frames and their products.
		Min_size refers to codon number instead of nucleotide number.
		Termination codons should map to an empty string."""
		if not isinstance(nested, bool):
			raise TypeError(f"The nested parameter is either True or False.\n")
		if not isinstance(min_size, int):
			raise TypeError(f"Codons are a quantum unit meaning that there are no half codons and none less that the quanta i.e no negative units in our case.\n")
		if not isinstance(start_codons, list):
			raise TypeError(f"The genetic code is a set of values, please use a list for the initiation codons.\n")
		if not isinstance(stop_codons, list):
			raise TypeError(f"The genetic code is a set of values, please use a list for the termination codons.\n")
		if not isinstance(codon_size, int):
			raise TypeError(f"The codon size should be the number of nucleotides per codon so the same criteria apply to the codon_size as teh min_size.\n NOTE: Please use the length of the codons used in the genetic_code dictionary for this variable.\n")
		if min_size < 0:
			raise ValueError(f"There is no negative open reading frame.\n")
		if codon_size != len(start_codons[0]):
			raise ValueError(f"The codon_size variable should be the same as the lenght of the codons used in the genetic code.\n")
		def holder(nested:bool = nested):
			orf = ""
			protein = ""
			self.orfs = []
			for reading_frame in range(codon_size):
				saving = False
				for _ in range(reading_frame, len(self.seq), codon_size):
					codon = self.seq[_:_ + codon_size:1] 
					if len(codon) != codon_size:
						return self.orfs
					if codon in start_codons: 
						start_codon_position = _
						saving=True
					if saving:
						orf += codon 
						protein+= genetic_code[codon]
					if codon in stop_codons:
						saving = False
						if nested:
							_ = start_codon_position
						if len(orf) >= codon_size * min_size:
							self.orfs.append([orf, protein])
						orf = ""
			return self.orfs
		if nested:
			return holder(nested = True)
		else:
			return holder(nested = False)

	def count_contents(self):
		"""Returns a list of the count of each nucleotide in the sequence."""
		counts = [["A", "T", "G", "C"], [0, 0, 0, 0]]
		for nuc in self.seq:
			for ref in range(4):
				if nuc == counts[0][ref]:
					counts[1][ref] += 1
		a = counts[1][0]
		t = counts[1][1]
		g = counts[1][2]
		c = counts[1][3]
		counts[1] = [a, t, g, c]
		return counts[1]
	
	def adenine_count(self):
		return self.count_contents()[0]
	
	def thymine_count(self):
		return self.count_contents()[1]
	
	def guanine_count(self):
		return self.count_contents()[2]
	
	def cytocine_count(self):
		return self.count_contents()[3]

	def percent_contents(self):
		"""Returns a list of the percentage of each nucleotide in the sequence."""
		counts = [["A", "T", "G", "C"], [0, 0, 0, 0]]
		for nuc in self.seq:
			for ref in range(4):
				if nuc == counts[0][ref]:
					counts[1][ref] += 1
		a = counts[1][0]
		t = counts[1][1]
		g = counts[1][2]
		c = counts[1][3]
		counts[1] = [a/len(self.seq)*100, t/len(self.seq)*100, g/len(self.seq)*100, c/len(self.seq)*100]
		return counts[1]

	def adenine_percent(self):
		return self.percent_contents()[0]
	
	def thymine_percent(self):
		return self.percent_contents()[1]
	
	def guanine_percent(self):
		return self.percent_contents()[2]
	
	def cytocine_percent(self):
		return self.percent_contents()[3]
	
class Gene(Sequence):
	"""This class is for representing genes, derived from the Sequence class."""
	def __init__(self, sequence: str|Sequence, exons: list = [[], []], name: str = "YES", description: str = "No description has been set", coding: bool = True):
		super().__init__(sequence)
		if not self.dna:
			raise ValueError(f"Please use DNA characters.\n")
		if not isinstance(name, str):
			raise TypeError(f"The name must be a string.\nExpecter type (str).\nGot{type(name)}.\n")
		if not isinstance(exons, list):
			raise TypeError(f"invalid exon type.\nExpected type (list).\nGot {type(exons)}")
		if not isinstance(description, str):
			raise TypeError(f"The description shiuld be a string.\nGiven type: {type(description)}")
		if not isinstance(coding, bool):
			raise TypeError(f"Please specify if this is the coding strand as True/False.\n")
		self.exons = exons
		self.introni()
		self.intron_seq()
		self.utr5 = self.exons[0][0]
		self.utr3 = self.exons[len(self.exons)-1][1]
		self.utrs = [self.seq[:self.exons[0][0]:1], self.seq[self.exons[len(self.exons)-1][1]::1]]
		self.coding = coding
		self.name = name
		self.lenght = len(self.seq)
		self.description = description
		self.nc_strand = self.reverse_complement()

	def __str__(self):
		"""Returns a format string containing all the gene's information."""
		return f"Sequence: {self.seq}, introns: {self.introns}, Introns: {self.introns}, Name: {self.name}, Length: {self.length}"
	
	def __len__(self):
		"""Returns a number representing the lenght of the gene."""
		return self.gene_length
		
	def transcribe(self):
		"""Returned the sequence's transcription product."""
		if self.coding:
			new = self.reverse_complement()
		else:
			new = self.reverse_complement()
		return new.replace("T", "U")

	def translate(self, code: dict = {"ATG": "M"}):
		"""Returns the product of the gene."""
		if not isinstance(code, dict):
			raise TypeError(f"The genetic code takes a dictionary of codons and amino-acids as input.\nExpected: type(dict).\nGot: {type(code)} instead.\n")
		product = ""

	def mature(self) -> list:
		"""Returns a matured version of the transcribed sequence.
		A matured version of the transcribed sequence is one where all the introns of the sequence have been excised."""
		new = [] 
		for val in self.exons:
			new.append(self.seq[val[0]:val[1]:1])
		return new

	def introni(self) -> list:
		"""Returns only the intronic indexes.
		Will not recognize twintrons."""
		self.introns = []
		for i in range(len(self.exons) - 1):
			try:
				self.introns.append(self.exons[i][1])
				self.introns.append(self.exons[i+1][0])
			except IndexError:
				break
		return self.introns

	def intron_seq(self) -> list:
		"""Will return the DNA or RNA of the parts of the sequence that correspond to introns within a list."""
		self.intron_sequences = []
		for intr in range(0, len(self.introns) - 1, 2):
			self.intron_sequences.append(self.seq[self.introns[intr]:self.introns[intr+1]:1])
		return self.intron_sequences

class Intergenic(Sequence):
	"""This class represents inter-genic sequences, derived from the Sequence class."""
	def __init__(self, sequence, name = ""):
		super().__init__(sequence)
		if not isinstance(name, str):
			raise ValueError(f"The name must be a string.\nExpecter type (str).\nGot{type(name)}.\n")
		self.name = name
		self.length = len(self.seq)

	def __str__(self):
		"""Returns a formated string of the sequence and name of the intra-genic sequence."""
		return f"Lenght: {self.length}, Name: {self.name}"
	
	def __len__(self):
		"""Returns the lenght of the sequence."""
		return self.lenght

class Chromosome(Sequence):
	"""This class represents a chromosome, chromosomes are sequences and thus this class inherits from the Sequence class."""
	def __init__(self, sequence: str|Sequence, genes: list, centromer: int|float, name = ""):
		super().__init__(sequence)
		self.genes = genes 
		self.intr = []
		self.chromosome_name = name
		self.center = centromer
	
	def __str__(self):
		"""Returns the sequence of the instance."""
		return self.seq
	
	def __len__(self):
		"""Returns the lenght of the chromosome's sequence."""
		return len(self.seq) 

class Genome(Sequence):
	"""This class represents a genome, genomes are made of sequences and if you combine all the sequences you get one sequence for the entires genome thus this class inherits from the sequence class."""
	def __init__(self, sequence, cpairs = [[[300, 400, 500], "1", 0.7, "AAT", True, False], []]) -> None:
		super().__init__(sequence)
		self.length = cpairs

	def __str__(self) -> int:
		"""Returns the sets of chromosomes."""
		return str(self.length)
	
	def __len__(self) -> int:
		"""Returns the genome lenght."""
		return self.length
		
class Person(Sequence):
	"""This class represents one single person."""
	def __init__(self, sequence = random.choices(population = ["A", "T", "G", "C"], cum_weights = (25, 50, 75, 100), k = len("AAAGGTACGCGCGCCGGCGCGTATAGCTTTAGTCGTGGACGCTAGCTAGCTGGTAGCGACAGGCGAGAAATGCTAGCATCGAGCATGCAGCGTTC")), genome = 10, ID = None, x_pos = 0, y_pos = 0) -> None:
			if not isinstance(genome, int|list):
				raise ValueError(f"The genome must be a list containing numbers.\nExpected something like {[1, 2, 3, 4, 5, 6, 7, 7, 8, 9, 10]} got {genome}.\n")
			if ID != None and not isinstance(ID, int):
				raise TypeError (f"The ID must be a natural number.\nExpected: type(int).\nGot: {type(ID)}.\n")
			if isinstance(ID, int):
				if 0 <= ID:
					pass
				else:
					raise TypeError(f"Improper value for ID.\nMust be a natural number.\n")
			if not isinstance(sequence, str|list):
				raise TypeError(f"The sequence must be a list of the nucleotides.\nExpected: type(list).\nGot: {type(sequence)}.\n")
			if not isinstance(x_pos, int) or x_pos < 0:
				raise TypeError(f"The positions of the individual start form the left most top postition of the window thus they cannot be negative and not real numbers.\n")
			if not isinstance(y_pos, int) or y_pos < 0:
				raise TypeError(f"The positions of the individual start form the left most top postition of the window thus they cannot be negative and not real numbers.\n")
			super().__init__(sequence)
			# Parametising the sequence
			self.genome = genome
			# Assigning a unique code to the individual
			self.id = ID
			# Coordinates
			self.x_pos = x_pos
			self.y_pos = y_pos
			# Setting status
			self.living = True
			# Set an image that represent the phenotype
			os.mkdir("./tmp/", exist_ok = true)
			self.phenotype = open(f"./tmp/person {self.id}.ppm", "w")
			self.phenotype.write(f"P3\n3 3\n")
			self.phenotype.close()

	def __str__(self) -> int:
		"""Returns the individuals sequence."""
		return self.seq

	def __len__(self) -> int:
		"""Returns the persons sequence lenght."""
		return len(self.seq)

	def mutate(self, rate: int|float = 0.5) -> None:
		"""Creates multiple point mutatuions for the person's sequnece using a rate instead of an ammount.
		Alters the original sequence."""
		self.mutations(ammount = math.ceil(len(self.seq) * rate), mutation_weights = {"point mutations": 20, "transversions": 40, "transitions": 60, "insertions": 60, "deletions": 60, "translocations": 80}) 

	def move(self, x_ammount: int, y_ammount: int) -> None:
		"""Changes the position of the individual."""
		if not isinstance(x_ammount, int) or x_ammount < - self.x_pos:
			raise TypeError(f"The x_pos cannot be negative or a non integer meaning that the x_ammount must not be a non integer and not less than the negative ammount of x_pos (so that x_pos + x_ammount !<0)")
		if not isinstance(y_ammount, int) or y_ammount < - self.y_pos:
			raise TypeError(f"The y_pos cannot be negative or a non integer meaning that the y_ammount must not be a non integer and not less than the negative ammount of y_pos (so that y_pos + y_ammount !<0)")
		self.x_pos += x_ammount
		self.y_pos += y_ammount 

	def make_window(self) -> None:
		"""Creates a window with a grid representation of the persons position."""
		window = Tk()
		window.title("R.S.R. Simulation")
		frame = ttk.Frame(window, padding = 20)
		frame.grid()
		window.mainloop()

class Population:
	def __init__(self, size = 100, generations = 100, reference_sequence = "AAAGGTACGCGCGCCGGCGCGTATAGCTTTAGTCGTGGACGCTAGCTAGCTGGTAGCGACAGGCGAGAAATGCTAGCATCGAGCATGCAGCGTTC", **instructions):
		if isinstance(reference_sequence, Sequence):
			self.reference = reference_sequence.seq	
		else:
			if not isinstance(reference_sequence, str|list):
				raise ValueError(f"Your reference sequence is invalid.\n")
			else:
				self.reference = Sequence(reference_sequence)
		if not isinstance(generations, int) or generations <= 0:
			raise ValueError(f"Please input a valid value for the generations.\nMust be an integer: Expected number > 0.\nGot {generations}.\n")
		if not isinstance(size, int):
			raise ValueError(f"The size must be a positive number.\nExpected int within range(1 to + infinity), instead got {size}.\n")
		# Number of people to start with
		self.size = size
		# Number of generations
		self.gen = generations
		# A list that keeps all the individuals (as the data structure above)
		self.people = []
		# list of all generations
		self.generations = []
		# String where we output all the data of the seqeunces
		self.log = f""
		# A list that keep track of all the generations individual count
		self.slog = [self.size]
		# A list of the distances from the reference
		self.distances = []
		# A list of the average reference similarity of a generation for each generation
		self.average = []
		# Tkinter stats
		self.window = Tk()
		self.window.geometry("1920x1080")
		self.window.title("R.S.R. Simulation")
		self.frame = ttk.Frame(self.window, padding = 10)
		self.frame.grid()
		ttk.Button(self.window, text = "Quit", command = self.window.destroy).grid(column = 10, row = 10)
		# A loop to create all the individuals we start with
		for _ in range(size):
			person = Person(sequence = random.choices(population = ["A", "T", "G", "C"], cum_weights = [self.reference.adenine_percent(), self.reference.adenine_percent() + self.reference.thymine_percent(), self.reference.adenine_percent() + self.reference.thymine_percent() + self.reference.guanine_percent(), self.reference.adenine_percent() + self.reference.thymine_percent() + self.reference.guanine_percent() + self.reference.cytocine_percent()], k = len(self.reference)), ID = _ + 1, x_pos = random.randint(0, 10), y_pos = random.randint(0, 10))
			# Shouldn't it be person.seq?
			self.people.append(person)
			rounded_dist = round(compare(self.reference.seq, person.seq), 2)
			self.log += f"Person {_ + 1}, Generation: Parent, {str(person)}, reference similarity ({self.reference}) roughly {rounded_dist}, No ancestors\n"
			self.distances.append(rounded_dist)
		self.average.append(round(stat.mean(self.distances), 2))
		# We need more variabbles to be able to pass the recursive structure into the bigtree functions
		# Family tree parameter
		self.ftree = ""
		self.make_window()

	def __str__(self) -> str:
		"""Returns a log of all the people that are and ever have been in the population."""
		return self.log

	def __len__(self) -> int:
		"""Returns the current ammount of individuals in the population."""
		return self.size

	def cross(self, verbose: bool = False, log: bool = False, criteria: int|float = 10, drift: int|float = 10) -> None:
		"""Picks two random individuals from the population and recombines their genome at a random location."""
		# cross refers more so to a recombination for DNA of two individuals
		if not isinstance(criteria, int|float):
			raise TypeError(f"Please give slection criteria.\nGot {criteria}.\n")
		if not isinstance(drift, int|float):
			raise TypeError(f"Please give a valid drift.\nGot {drift}.\n")
		if not isinstance(verbose, bool):
			raise TypeError(f"Problem with given argument.\nVerbos must be boolean, got {verbose}.\n")
		if not isinstance(log, bool):
			raise TypeError(f"Invalid argument for log. Given {log}.\n")
		if drift < 0:
			raise ValueError(f"Drift does not take negative values.\n")
		# Why do we need this?
		self.distances = []
		self.ancestors = []
		new = []
		toremovelist = [] 
		for loop in range(self.size * math.ceil(abs((self.slog[0] - self.size)/self.slog[len(self.slog) - 2])) + 5): 
			person1 = random.choice(self.people)
			person2 = random.choice(self.people)
			while person1 == person2 and self.size > 1:
				person2 = random.choice(self.people)
			rec_position = random.randint(0, len(str(person1)) - 1)
			person3 = Person(sequence = str(person1)[:rec_position:1] + str(person2)[rec_position::1], ID = loop + 1, x_pos = math.ceil((person2.x_pos + person1.x_pos)/2), y_pos = math.ceil((person2.y_pos + person1.y_pos)/2))
			if verbose:
				print(f"Person {person1.id}: {str(person1)} x Person {person2.id}: {str(person2)} at position {rec_position} To give new Person {person3.id}: {str(person3)}.")
			if log:
				self.log += (f"Person {person1.id}: {str(person1)} x Person {person2.id}: {str(person2)} at position {rec_position} To give new Person {person3.id}: {str(person3)}.\n")
			# Now self.ancestors becomes a list with people when before it had sequences
			self.ancestors.append(person1)
			self.ancestors.append(person2)
			# New containt a person but then self.sequences turns to new so seqeunces does not contain sequences
			new.append(person3)
			self.distances.append(round(compare(self.reference.seq, person3.seq), 2))
		if log:
			self.log += f"\n"
		self.sequences = new
		self.size = len(self.sequences)
		self.average.append(round(stat.mean(self.distances), 2))
		self.history.append(self.ancestors)
		for person in self.sequences:
			if random.randint(0, 100) < drift or compare(person.seq, self.reference.seq) <= self.average[len(self.average) - 2] - criteria:
				if verbose:
					print(f"Person {person.id} {person} has not reached adulthood.")
				if log:
					self.log += (f"Person {person.id} {person} has not reached adulthood.\n")
				toremovelist.append(person)
			else:
				if verbose:
					print(f"Person {person.id} {person.seq} has reached adulthood.")
				if log:
					self.log += f"Person {person.id} {person.seq} has reached adulthood.\n"
		for ind in toremovelist:
			"add the removal from the gui in here"
			self.sequences.remove(ind)
		self.size = len(self.sequences)

	def mutate(self, rate: int|float = 0.5, verbose: bool = False, log: bool = False, criteria: int|float = 10, drift: int|float = 10) -> None:
		"""Mutates every individual of the population."""
		if not isinstance(rate, int|float):
			raise TypeError(f"Please give numeral rate.\nGot {rate}.\n")
		if not isinstance(criteria, int|float):
			raise TypeError(f"Please give slection criteria.\nGot {criteria}.\n")
		if not isinstance(drift, int|float):
			raise TypeError(f"Please give a valid drift.\nGot {drift}.\n")
		if not isinstance(verbose, bool):
			raise TypeError(f"Problem with given argument.\nVerbose must be boolean, got {verbose}.\n")
		if not isinstance(log, bool):
			raise TypeError(f"Invalid argument for log. Given {log}.\n")
		if drift < 0:
			raise ValueError(f"Drift does not take negative values.\n")
		if rate <= 0:
			raise ValueError(f"rate must be positive non 0.\nGiven {rate}.\n")
		self.distances = []
		for person in self.sequences:
			prevseq = person.seq
			person.mutate(rate = rate)
			self.distances.append(round(compare(self.reference.seq, person.seq), 2))
			if verbose:
				print(f"Person {person.id} {prevseq} has mutated to {person}.")
			if log:
				self.log += f"Person {person.id} {prevseq} has mutated to {person}.\n"
		if log:
			self.log += f"\n"                                               
		# This seems redundant why cant i d one loop for all the selection? If not this way the loop number is variable and the results are un predictable 
		# However if some smarter method for looping the individuals where to be made so that the loop can continue with the changing number
		toremovelist = [] 
		for person in self.sequences:
			if random.randint(0, 100) < drift or compare(person.seq, self.reference.seq) <= self.average[len(self.average) - 2] - criteria:
				if verbose:
					print(f"Person {person.id} {person} has not lived long enough to see the next generation.")
				if log:
					self.log += (f"Person {person.id} {person} has not lived long enough to see the next generation.\n")
				toremovelist.append(person)
			else:
				if verbose:
					print(f"Person {person.id} {person.seq} has lived long enough to see the next generation.")
				if log:
					self.log += f"Person {person.id} {person.seq} has lived long enough to see the next generation.\n"
		for ind in toremovelist:
			self.sequences.remove(ind)
		self.size = len(self.sequences)

	def family_tree(self,verbose: bool = False, log: bool = True):
		if not isinstance(log, bool):
			raise TypeError(f"Problem with given argument.\nLog must be boolean, got {log}.\n")
		if not isinstance(verbose, bool):
			raise TypeError(f"Problem with given argument.\nVerbose must be boolean, got {verbose}.\n")
		# self.ancestors worn't work here
		# self.ftree = bt.list_to_tree(self.ancestors)
		return self.ftree

	def generations(self, rate: int|float = 0.5, drift: int|float = 10, criteria: int|float = 10, verbose: bool = False, log: bool = True) -> tuple:
		"""Simulates te entire population throughout all of its life span."""
		if not isinstance(drift, int|float):
			raise TypeError(f"Please give a valid drift.\nGot {drift}.\n")
		if not isinstance(criteria, int|float):
			raise TypeError(f"Please give slection criteria.\nGot {criteria}.\n")
		if not isinstance(rate, int|float):
			raise TypeError(f"Please give numeral rate.\nGot {rate}.\n")
		if not isinstance(verbose, bool):
			raise TypeError(f"Problem with given argument.\nVerbose must be boolean, got {verbose}.\n")
		if not isinstance(log, bool):
			raise TypeError(f"Problem with given argument.\nLog must be boolean, got {log}.\n")
		for cgen in range(self.gen):
			# This is nice in the output but ugly in the code
			info = f"""
# ######################################################################

# Generation: {cgen + 1}

# Population size: {self.size}

# Average population reference similarity: {self.average[cgen]}

# ######################################################################\n"""
			if verbose:
				print(info)
			if log:
				self.log += info + "\n" 
			self.mutate(rate = rate, verbose = verbose, log = log, criteria = criteria, drift = drift)
			for widget in self.frame.winfo_children():
				widget.destroy()
			self.make_window()
			# Multiple if verbose
			if verbose:
				print(self.family_tree(verbose = verbose, log = log))
			# Two if self.size <= 1?
			# Both cross and mutate can kill the population hence the double if
			if self.size <= 1:
				self.slog.append(self.size)
				print("The populations has died.\n")
				self.log += "The populations has died."
				return self.slog, self.average, self.log
			self.cross(population = self.sequences, verbose = verbose, log = log, criteria = criteria, drift = drift)
			self.slog.append(self.size)
			if self.size <= 1:
				print("The populations has died.\n")
				self.log += "The populations has died."
				return self.slog, self.average, self.log
		return self.slog, self.average, self.log
		
	def make_window(self) -> None:
		"""Simulates the populations movements."""
		for person in self.sequences:
			person.move(x_ammount = random.randint(0, 10), y_ammount = random.randint(0, 10))
			ttk.Label(self.frame, text = person.id).grid(column = person.x_pos, row = person.y_pos)
		self.window.update()

def cmd_line_input():
	"""Returns all command line argument inputs as variables for the program to use."""
	# There must be a better way to add arguments
	try:
		parameters = {}
		for i in range(1, len(sys.argv)):
			arg = sys.argv[i]
			if arg == "-o" or arg == "-acc":
				parameters[arg] = sys.argv[i + 1]
			elif arg == "-size" or arg == "-gen":
				parameters[arg] = int(sys.argv[i + 1])
			elif arg == "-rate":
				parameters[arg] = float(sys.argv[i + 1])
			elif arg == "-h":
				sys.exit("Usage:\n[-h][shows a list of all commands].\n[-acc][accession number: can be empty].\n[-email][Your email. Required if an accession number is provided to fetch the sequence from NCBI].\n[-o][output_path: can be empty (inputs default value)].\n[-gen][number for size of each generation: Must be greater than 0. Can be empty (inputs default value)].\n[-size][number for size of the population: Must be greater than 0. Can be empty (inputs default value)].\n[-rate][mutation rate: Must be greater than 0. Can be empty (inputs default value)].\n[-v][Verbose: prints processes on screen].\n[-q][Quiet: stops printing output (may increase speed)].\n[-drift][Defines a genetic drift parameter that kills individuals randomly].\n[-criteria][Defines a selection criteria for natural selection, warning may cause population to die easily].\n[-log][Creates a log file to track run].\n")
				#sys.exit("Usage:\n[-h][shows a list of all commands].\n[-acc][accession number: can be empty].\n[-o][output_path: can be empty (inputs default value)].\n[-gen][number for size of each generation: Must be greater than 0. Can be empty (inputs default value)].\n[-size][number for size of the population: Must be greater than 0. Can be empty (inputs default value)].\n[-rate][mutation rate: Must be greater than 0. Can be empty (inputs default value)].\n[-v][Verbose: prints processes on screen].\n[-q][Quiet: stops printing output (may increase speed)].\n[-speed][Speed at which the program runs: An exponential scaling speed option, note that any value lower than the default may cause infinite run time. Values range from 1 to 100].\n[-regression][Allows the reference similarity to drop during the run].\n[-drift][Defines a genetic drift parameter that kills individuals randomly].\n[-criteria][Defines a selection criteria for natural selection, warning may cause population to die easily].\n[-log][Creates a log file to track run].\n")
			elif arg == "-v":
				parameters[arg] = True
				parameters["-q"] = False
			elif arg == "-q":
				parameters[arg] = True
				parameters["-v"] = False
			# elif arg == "-speed":
			# 	parameters[arg] = int(sys.argv[i + 1])
			# elif arg == "-regression":
			# 	parameters[arg] = True
			elif arg == "-log":
				parameters[arg] = True
			elif arg == "-drift":
				parameters[arg] = int(sys.argv[i + 1])
			elif arg == "-criteria":
				parameters[arg] = int(sys.argv[i + 1])
			elif arg == "-email":
				parameters[arg] = str(sys.argv[i + 1])
		# setting default values
		if "-rate" not in sys.argv:
			parameters["-rate"] = 0.5
		if "-email" not in sys.argv:
			parameters["-email"] = "a.Podotas@gmail.com"
		if "-o" not in sys.argv:
			parameters["-o"] = "."
		if "-size" not in sys.argv:
			parameters["-size"] = 100
		if "-gen" not in sys.argv:
			parameters["-gen"] = 100
		if "-acc" not in sys.argv:
			parameters["-acc"] = "AAAGGTACGCGCGCCGGCGCGTATAGCTTTAGTCGTGGACGCTAGCTAGCTGGTAGCGACAGGCGAGAAATGCTAGCATCGAGCATGCAGCGTTC"
		if "-q" not in sys.argv and "-v" not in sys.argv:
			parameters["-v"] = False
			parameters["-q"] = True
		elif "-q" in sys.argv and "-v" not in sys.argv:
			parameters["-v"] = False
			parameters["-q"] = True
		elif "-v" in sys.argv and "-q" in sys.argv:
			sys.exit(f"Can't have -q and -v together.\n")
		else:
			parameters["-v"] = True
			parameters["-q"] = False
		# if "-regression" not in sys.argv:
		# 	parameters["-regression"] = False
		# if "-speed" not in sys.argv:
		# 	parameters["-speed"] = 100
		if "-log" not in sys.argv:
			parameters["-log"] = True
		if "-drift" not in sys.argv:
			parameters["-drift"] = 10
		if "-criteria" not in sys.argv:
			parameters["-criteria"] = 10
	except(ValueError):
		sys.exit(f"Problem with a parameter value.\nUsage:\n[-h][shows a list of all commands].\n[-acc][accession number: can be empty].\n[-email][Your email. Required if an accession number is provided to fetch the sequence from NCBI].\n[-o][output_path: can be empty (inputs default value)].\n[-gen][number of generations: Must be greater than 0. Can be empty (inputs default value)].\n[-size][number for size of each generation: Must be greater than 0. Can be empty (inputs default value)].\n[-rate][mutation rate: Must be greater than 0. Can be empty (inputs default value)].\n[-v][Verbose: prints processes on screen].\n[-q][Quiet: stops printing output (may increase speed)].\n[-drift][Defines a genetic drift parameter that kills individuals randomly].\n[-criteria][Defines a selection criteria for natural selection, warning may cause population to die easily].\n[-log][Creates a log file to track run].\nProblem in parameter value.\nGiven{sys.argv}.\n")
		#sys.exit(f"Problem with a parameter value.\nUsage:\n[-h][shows a list of all commands].\n[-acc][accession number: can be empty].\n[-o][output_path: can be empty (inputs default value)].\n[-gen][number of generations: Must be greater than 0. Can be empty (inputs default value)].\n[-size][number for size of each generation: Must be greater than 0. Can be empty (inputs default value)].\n[-rate][mutation rate: Must be greater than 0. Can be empty (inputs default value)].\n[-v][Verbose: prints processes on screen].\n[-q][Quiet: stops printing output (may increase speed)].\n[-speed][Speed at which the program runs: An exponential scaling speed option, note that any value lower than the default may cause infinite run time. Values range from 1 to 100].\n[-regression][Allows the reference similarity to drop during the run].\n[-drift][Defines a genetic drift parameter that kills individuals randomly].\n[-criteria][Defines a selection criteria for natural selection, warning may cause population to die easily].\n[-log][Creates a log file to track run].\nProblem in parameter value.\nGiven{sys.argv}.\n")
	for key in sys.argv[1::1]:
		if key in parameters.keys() or key in str(parameters.values()) or key == ".\\test_project.py" or key == "test_project.py": 
			continue
		else:
			sys.exit(f"Promblem with given argument.\nUsage:\n[-h][shows a list of all commands].\n[-acc][accession number: can be empty].\n[-email][Your email. Required if an accession number is provided to fetch the sequence from NCBI].\n[-o][output_path: can be empty (inputs default value)].\n[-gen][number of generations: Must be greater than 0. Can be empty (inputs default value)].\n[-size][number for size of each generation: Must be greater than 0. Can be empty (inputs default value)].\n[-rate][mutation rate: Must be greater than 0. Can be empty (inputs default value)].\n[-v][Verbose: prints processes on screen].\n[-q][Quiet: stops printing output (may increase speed)].\n[-drift][Defines a genetic drift parameter that kills individuals randomly].\n[-criteria][Defines a selection criteria for natural selection, warning may cause population to die easily].\n[-log][Creates a log file to track run].\nProblem in arguments.\nGiven{sys.argv}.\nUnable to interpret [{key}].\n")
			#sys.exit(f"Promblem with given argument.\nUsage:\n[-h][shows a list of all commands].\n[-acc][accession number: can be empty].\n[-o][output_path: can be empty (inputs default value)].\n[-gen][number of generations: Must be greater than 0. Can be empty (inputs default value)].\n[-size][number for size of each generation: Must be greater than 0. Can be empty (inputs default value)].\n[-rate][mutation rate: Must be greater than 0. Can be empty (inputs default value)].\n[-v][Verbose: prints processes on screen].\n[-q][Quiet: stops printing output (may increase speed)].\n[-speed][Speed at which the program runs: An exponential scaling speed option, note that any value lower than the default may cause infinite run time. Values range from 1 to 100].\n[-regression][Allows the reference similarity to drop during the run].\n[-drift][Defines a genetic drift parameter that kills individuals randomly].\n[-criteria][Defines a selection criteria for natural selection, warning may cause population to die easily].\n[-log][Creates a log file to track run].\nProblem in arguments.\nGiven{sys.argv}.\nUnable to interpret [{key}].\n")
	return parameters

def NCBI_parse(accession: str, email: str):
	"""Returns the sequence of the command line as a usable variable and interprits NCBI accession sequences."""
	# in case the input is a literal sequence
	if not isinstance(accession, str):
		raise TypeError(f"Please give a valid accession.\nGot {accession}.\n")
	if re.search(r"[AGCT]+", accession):
		return str(accession)
	# in case of -acc being an NCBI accession number, fetching sequence
	elif re.search(r"(\w\d{5}.{0,3})|(\w{2}\d{6}.{0,3})", accession):
		Entrez.email = email
		entry = Entrez.efetch(db = "nucleotide", id = accession, retmode = "text", rettype = "gb")
		for query in SeqIO.parse(entry, "gb"):
			seq = str(query.seq)
		return seq
	else:
		raise ValueError(f"Either a sequence or accession number must be input\nGiven: {accession}")

def make_folders(output_path: str, folder_name: str):
	"""Creates the folder for all the ouput files as specified by the command line."""
	if not isinstance(output_path, str):
		raise TypeError(f"Please specify path.\nGot {output_path}.\n")
	if not isinstance(folder_name, str):
		raise TypeError(f"Please specify path.\nGot {folder_name}.\n")
	print(f"""
##########################################################################

Creating output folder

##########################################################################\n""")
	new_folder_path = os.path.join(output_path, folder_name)
	os.makedirs(new_folder_path, exist_ok = True)

# Make it so that it uses the needleman wunsch
def compare(seq1: Sequence|str|list|tuple|set, seq2: Sequence|str|list|tuple|set):
	"""Sequence aligment based on Needleman - Wunsch"""
	# calculates the percent identity of two sequences
	if not isinstance(seq2, Sequence|str|list|tuple|set):
		raise TypeError(f"Cannot calculate distance for {seq2}.\n")
	if not isinstance(seq1, Sequence|str|list|tuple|set):
		raise TypeError(f"Cannot calculate distance for {seq1}.\n")
	if len(seq1) == len(seq2):
		matches = 0
		for nucleotide in range(len(seq1)):
			if seq1[nucleotide] == seq2[nucleotide]:
				matches += 1
		return (matches/len(seq1))*100
	# Make the cases exist
	elif len(seq1) > len(seq2):
		return 0
	else:
		return 0

def make_figures(size: list|tuple|set, ref: list|tuple|set):
	"""Creates a figure png image in the output folder representing what the generations function has done in its last run."""
	if not isinstance(size, list|tuple|set):
		raise TypeError(f"figures need iterable data for {size}.\n")
	if not isinstance(ref, list|tuple|set):
		raise TypeError(f"figures need iterable data for {ref}.\n")
	sz = plt.subplot(2, 1, 1)
	sz.plot(range(len(size)), size, label = "Population Size", c = "hotpink", ls = "-", marker = "o", alpha = 0.5)
	sz.set_xlim(-1, len(size))
	sz.set_ylim(0, max(size) + math.ceil(math.sqrt(max(size))))
	sz.set_ylabel("Population Size")
	sz.grid(ls = ":")
	sz.legend()
	refd = plt.subplot(2, 1, 2)
	refd.plot(range(len(ref)), ref, label = "Average Reference Similarity", c = '#4CAF50', ls = "-", marker = "o", alpha = 0.5)
	refd.set_xlim(-1, len(ref))
	refd.set_ylim(0, 110)
	refd.set_ylabel("Percent Identity")
	refd.set_xlabel("Generation")
	refd.grid(ls = ":")
	refd.legend()
	plt.suptitle("RSR Output")
	plt.savefig(f"{cmd_line_input()['-o']}/RSR_output_folder/RSR Output Graphs.png")
	return True

def main():
	main_time = time.time()
	cmd = cmd_line_input()
	make_folders(cmd['-o'], "RSR_output_folder")
	athens = Population(size = cmd["-size"], generations = cmd["-gen"], reference_sequence = NCBI_parse(cmd["-acc"], cmd["-email"]))
	output_file = os.path.join("RSR_output_folder", 'Results.txt')
	generations_output = athens.generations(rate = cmd["-rate"], verbose = cmd["-v"], log = cmd["-log"], drift = cmd["-drift"], criteria = cmd["-criteria"])
	with open(output_file, 'w') as out_file:
		out_file.write(f"Reference sequence: {NCBI_parse(cmd['-acc'], cmd['-email'])}\n\nPopulation:\n{generations_output[2]}")
		make_figures(generations_output[0], generations_output[1])
	print(f"""
##############################################################################################################################################

Process finished.\n\nParse results at {os.path.abspath(cmd['-o'])}\RSR_output_folder within Results.txt and RSR Outpout Graph.png

##############################################################################################################################################\n""")
	main_end_time = time.time() - main_time
	print(f"Took {main_end_time:.2f} seconds to run.\n")

if __name__ == "__main__":
	main()

