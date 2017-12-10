import simulation as sim
import os
import numpy as np

INI_file = "params.ini" #sys.argv[1]
output_dir = "output" #sys.argv[2]

def fitness (profile, target_profile):
    #difference entre profil et profil cible
    dif = profile/np.sum(profile) - target_profile
    # fitness = 1/distance
    return 1/np.sum(np.multiply(dif,dif))

def read_positions () :
    #lire les fichiers de tousgenesidentiques
    config = sim.read_config_file(INI_file)
    TSS_file = config.get('INPUTS', 'TSS')
    TTS_file = config.get('INPUTS', 'TTS')
    tss = sim.load_gff(TSS_file)
    tts = sim.load_tab_file(TTS_file)
    pos_start = tss['TSS_pos'].values
    pos_end = tts['TTS_pos'].values
    orientation  = sim.str2num(tss['TUorient'].values)
    return [[pos_start[i],pos_end[i], orientation[i]] for i in range(len(pos_start))]

def write_positions () :
    #écrire les fichiers de tousgenesidentiques (TSS, TTS et GFF)
    return []

def mutations (genome) :
    positions = list(genome)
    #inversion

    #deletion

    #insertion

    write_positions(positions)
    return positions

def metropolis() :
    envir_file = open(os.path.join('tousgenesidentiques', 'environment.dat'))
    target_profile = []
    for line in envir_file.readlines():
        target_profile.append(float(line.split()[1]))
    target_profile = np.array(target_profile)
    envir_file.close()
    #simulation :
    current_genome = read_positions()
    transcriptome = np.array(sim.start_transcribing(INI_file, output_dir))
    current_fitness = fitness(transcriptome, target_profile)
    hist_fitness = [current_fitness]

    for i in range(2) :
        current_genome = mutations(current_genome)
        transcriptome = np.array(sim.start_transcribing(INI_file,output_dir))
        current_fitness = fitness(transcriptome,target_profile)
        print(current_fitness, transcriptome)
        hist_fitness.append(current_fitness)

    print(current_fitness)
    return []

metropolis()