import Evo_RSR as rsr
import time

#Alternatively import the moduels from requierments individually with Evo_RSR

def main():
	main_time = time.time()
	cmd = rsr.cmd_line_input()
	rsr.make_folders(cmd['-o'], "RSR_output_folder")
	athens = rsr.Population(size = cmd["-size"], generations = cmd["-gen"], reference_sequence = rsr.NCBI_parse(cmd["-acc"]))
	output_file = rsr.os.path.join("RSR_output_folder", 'Results.txt')
	generations_output = athens.generations(rate = cmd["-rate"], verbose = cmd["-v"], log = cmd["-log"], drift = cmd["-drift"], criteria = cmd["-criteria"])
	with open(output_file, 'w') as out_file:
		out_file.write(f"Reference sequence: {rsr.NCBI_parse(cmd['-acc'])}\n\nPopulation:\n{generations_output[2]}")
		rsr.make_figures(generations_output[0], generations_output[1])
	print(f"""
##############################################################################################################################################

Process finished.\n\nParse results at {rsr.os.path.abspath(cmd['-o'])}\RSR_output_folder within Results.txt and RSR Outpout Graph.png

##############################################################################################################################################\n""")
	main_end_time = time.time() - main_time
	print(f"Took {main_end_time:.2f} seconds to run.\n")

if __name__ == "__main__":
	main()

