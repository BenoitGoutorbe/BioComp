import simulation as sim
import os
import numpy as np
import random

INI_file = "params.ini" #sys.argv[1]
output_dir = "output" #sys.argv[2]

config = sim.read_config_file(INI_file)
TSS_file = config.get('INPUTS','TSS')
TTS_file = config.get('INPUTS','TTS')
GFF_file = config.get('INPUTS', 'GFF')
P_DEL = 0.01
P_INS = 0.01
P_INV = 0.01

def fitness (profile, target_profile):
    #difference entre profil et profil cible
    dif = profile/np.sum(profile) - target_profile
    # fitness = 1/distance
    return 1/np.sum(np.multiply(dif,dif))

def read_positions () :
    #lire les fichiers de tousgenesidentiques
    tss = sim.load_gff(TSS_file)
    tts = sim.load_tab_file(TTS_file)
    pos_start = tss['TSS_pos'].values
    pos_end = tts['TTS_pos'].values
    orientation  = sim.str2num(tss['TUorient'].values)
    return [[pos_start[i],pos_end[i], orientation[i]] for i in range(len(pos_start))]

def write_positions (positions) :
    #écrire les fichiers de tousgenesidentiques (TSS, TTS et GFF)
    return []

def mutations (genome) :
    positions = list(genome)
    gff_df_raw = sim.load_gff(GFF_file)
    size_genome = int(list(gff_df_raw)[4])
    codant = []
    for g in positions :
        if g[2] > 0 :
            codant = codant + list(range(g[0],g[1]+1))
        else :
            codant = codant + list(range(g[0],g[1]-1, -1))
    non_codant= [i for i in range(1,size_genome+1) if i not in codant]
    print(len(codant), len(non_codant), non_codant[0],non_codant[-1])

    #inversion
    if random.random() < P_INS :
        

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