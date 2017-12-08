import sys
import simulation as sim

INI_file = "params.ini" #sys.argv[1]
output_dir = "output" #sys.argv[2]

sim.start_transcribing(INI_file, output_dir)
