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
P_DEL = 0.0
P_INS = 0.0
P_INV = 1.0
SIZE_INDEL = 1000


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

    #inversion
    if random.random() < P_INV :
        barrieres = random.sample(non_codant,2)
        barrieres.sort()
        print('INVERSION : intervalle', barrieres)
        for g in positions :
            if g[0] > barrieres[0] and g[1] < barrieres[1]:
                g[0] = barrieres[1] - (g[0] - barrieres[0])
                g[1] = barrieres[1] - (g[1] - barrieres[0])
                g[2] = -1.0 * g[2]

    #deletion
    if random.random() < P_DEL:
        position_del = random.choice(non_codant)
        while position_del+SIZE_INDEL in codant :
            position_del = random.choice(non_codant)
        print('DELETION', position_del)
        size_genome = size_genome - SIZE_INDEL
        for g in positions:
            if g[0] > position_del and g[1] > position_del:
                g[0] = g[0] - SIZE_INDEL
                g[1] = g[1] - SIZE_INDEL

    #insertion
    if random.random() < P_INS:
        position_ins = random.choice(non_codant)
        print('INSERTION', position_ins)
        size_genome = size_genome + SIZE_INDEL
        for g in positions:
            if g[0] > position_ins and g[1] > position_ins:
                g[0] = g[0] + SIZE_INDEL
                g[1] = g[1] + SIZE_INDEL
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
    print(current_fitness,transcriptome)
    print(current_genome)
    for i in range(2) :
        current_genome = mutations(current_genome)
        transcriptome = np.array(sim.start_transcribing(INI_file,output_dir))
        current_fitness = fitness(transcriptome,target_profile)
        print(current_fitness, transcriptome)
        print(current_genome)
        hist_fitness.append(current_fitness)
    return []

metropolis()