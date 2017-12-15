import simulation as sim
import os
import numpy as np
import random

INI_file = "params.ini" #sys.argv[1]
output_dir = "output" #sys.argv[2]
trial_dir = "essai" #sys.argv[2]


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

def write_positions (positions, size_genome, pos_bar) :
    #écrire les fichiers de tousgenesidentiques (TSS, TTS et GFF)
    write_TTS(positions)
    write_TSS(positions)
    write_tousgenesidentiques(positions, size_genome)
    write_prot(posbar)
    return []

def write_TTS(positions) :
    file = open(os.path.join(trial_dir,'TTS.dat'), 'w+')
    file.write("TUindex\tTUorient\tTTS_pos\tTTS_proba_off\n")
    for i, v in enumerate(positions):
        file.write(str(i)+ "\t")
        if(v[2] == +1):
            file.write("+\t")
        else:
            file.write("-\t")
        file.write(str(v[1]) + "\t" + "1.\n")
    file.close()
    return []
        
def write_TSS(positions) :
    file = open(os.path.join(trial_dir,'TSS.dat'), 'w+')
    file.write("TUindex\tTUorient\tTSS_pos\tTSS_strength\n")
    for i, v in enumerate(positions):
        file.write(str(i)+ "\t")
        if(v[2] == +1):
            file.write("+\t")
        else:
            file.write("-\t")
        file.write(str(v[0]) + "\t" + ".2\n")
    file.close()
    return []
        
def write_tousgenesidentiques(positions, size_genome) :
    file = open(os.path.join(trial_dir,'tousgenesidentiques.gff'), 'w+')
    file.write("##gff-version 3\n#!gff-spec-version 1.20\n#!processor NCBI annotwriter\n##sequence-region tousgenesidentiques 1 ")
    file.write(str(size_genome) + "\n")
    file.write("tousgenesidentiques\tRefSeq\tregion\t1\t"+ str(size_genome) + "\t.\t+\t.\tID=id0;Name=tousgenesidentiques\n")
    for i, v in enumerate(positions):
        file.write("tousgenesidentiques\tRefSeq\tgene\t" + str(v[0]) + "\t" + str(v[1]) + "\t.\t")
        if(v[2] == +1):
            file.write("+\t.\tID=g1;Name=g"+str(i+1)+"\n")
        else:
            file.write("-\t.\tID=g1;Name=g"+str(i+1)+"\n")
    file.close()
    return []

    
def write_prot(posbar) :
    file = open(os.path.join(trial_dir,'prot.dat'), 'w+')
    file.write("prot_name\tprot_pos\n")
    for i,v in enumerate(posbar):
        file.write("hns\t" + str(v) + "\n")
    file.close()
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

    write_positions(positions, size_genome)
    return positions, size_genome

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



write_positions(positions, size_genome)
