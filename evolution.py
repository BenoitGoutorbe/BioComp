import simulation as sim
import os
import numpy as np
import random

INI_file = "params.ini" #sys.argv[1]
output_dir = "output" #sys.argv[2]
trial_dir = "tousgenesidentiques" #sys.argv[2]


config = sim.read_config_file(INI_file)
TSS_file = config.get('INPUTS','TSS')
TTS_file = config.get('INPUTS','TTS')
GFF_file = config.get('INPUTS', 'GFF')
Prot_file = config.get('INPUTS', 'BARR_FIX')
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
    orientation = sim.str2num(tss['TUorient'].values)
    positions_genes = [[pos_start[i],pos_end[i],orientation[i]] for i in range(len(pos_start))]
    gff_df_raw = sim.load_gff(GFF_file)
    size_genome = int(list(gff_df_raw)[4])
    prot = sim.load_tab_file(Prot_file)
    positions_barrieres = prot['prot_pos'].values
    return (positions_genes, positions_barrieres, size_genome)


def write_positions (positions_genes, positions_barrieres, size_genome) :
    #écrire les fichiers de tousgenesidentiques (TSS, TTS et GFF)

    #TTS.dat
    file = open(os.path.join(trial_dir,'TTS.dat'), 'w+')
    file.write("TUindex\tTUorient\tTTS_pos\tTTS_proba_off\n")
    for i, v in enumerate(positions_genes):
        file.write(str(i)+ "\t")
        if(v[2] == +1):
            file.write("+\t")
        else:
            file.write("-\t")
        file.write(str(v[1]) + "\t" + "1.\n")
    file.close()
        
    #TSS.dat
    file = open(os.path.join(trial_dir,'TSS.dat'), 'w+')
    file.write("TUindex\tTUorient\tTSS_pos\tTSS_strength\n")
    for i, v in enumerate(positions_genes):
        file.write(str(i)+ "\t")
        if(v[2] == +1):
            file.write("+\t")
        else:
            file.write("-\t")
        file.write(str(v[0]) + "\t" +".2\n")
    file.close()
        
    #tousgenesidentiques.gff
    file = open(os.path.join(trial_dir,'tousgenesidentiques.gff'), 'w+')
    file.write("##gff-version 3\n#!gff-spec-version 1.20\n#!processor NCBI annotwriter\n##sequence-region tousgenesidentiques 1 ")
    file.write(str(size_genome) + "\n")
    file.write("tousgenesidentiques\tRefSeq\tregion\t1\t"+ str(size_genome) + "\t.\t+\t.\tID=id0;Name=tousgenesidentiques\n")
    for i, v in enumerate(positions_genes):
        file.write("tousgenesidentiques\tRefSeq\tgene\t" + str(v[0]) + "\t" + str(v[1]) + "\t.\t")
        if(v[2] == +1):
            file.write("+\t.\tID=g1;Name=g"+str(i+1)+"\n")
        else:
            file.write("-\t.\tID=g1;Name=g"+str(i+1)+"\n")
    file.close()

    return []

def mutations (genome,barrieres_topo, size_genome) :
    positions_genes = list(genome)
    positions_barrieres = list(barrieres_topo)
    codant = []
    for g in positions_genes :
        if g[2] > 0 :
            codant = codant + list(range(g[0],g[1]+1))
        else :
            codant = codant + list(range(g[0],g[1]-1, -1))
    non_codant= [i for i in range(1,size_genome+1) if i not in codant]

    #inversion
    if random.random() < P_INV :
        bornes = random.sample(non_codant,2)
        bornes.sort()
        print('INVERSION : intervalle', bornes)
        for g in positions_genes :
            if g[0] > bornes[0] and g[1] < bornes[1]:
                g[0] = bornes[1] - (g[0] - bornes[0])
                g[1] = bornes[1] - (g[1] - bornes[0])
                g[2] = -1.0 * g[2]
        for bar in positions_barrieres :
            if bar > bornes[0] and bar < bornes[1]:
                bar = bornes[1] - (bar - bornes[0])

    #deletion
    if random.random() < P_DEL:
        position_del = random.choice(non_codant)
        while position_del+SIZE_INDEL in codant :
            position_del = random.choice(non_codant)
        print('DELETION', position_del)
        size_genome = size_genome - SIZE_INDEL
        for g in positions_genes:
            if g[0] > position_del and g[1] > position_del:
                g[0] = g[0] - SIZE_INDEL
                g[1] = g[1] - SIZE_INDEL
        for bar in positions_barrieres:
            if bar > position_del :
                bar = bar - SIZE_INDEL
    #insertion
    if random.random() < P_INS:
        position_ins = random.choice(non_codant)
        print('INSERTION', position_ins)
        size_genome = size_genome + SIZE_INDEL
        for g in positions_genes:
            if g[0] > position_ins and g[1] > position_ins:
                g[0] = g[0] + SIZE_INDEL
                g[1] = g[1] + SIZE_INDEL
        for bar in positions_barrieres:
            if bar > position_ins:
                bar = bar + SIZE_INDEL
    #write_positions(positions_genes, positions_barrieres, size_genome)
    return positions_genes, positions_barrieres, size_genome


def metropolis() :
    envir_file = open(os.path.join('tousgenesidentiques', 'environment.dat'))
    target_profile = []
    for line in envir_file.readlines():
        target_profile.append(float(line.split()[1]))
    target_profile = np.array(target_profile)
    envir_file.close()
    #simulation :
    current_genome, current_barrieres, current_size = read_positions()
    transcriptome = np.array(sim.start_transcribing(INI_file, output_dir))
    current_fitness = fitness(transcriptome, target_profile)
    hist_fitness = [current_fitness]
    print(current_fitness,transcriptome)
    print(current_genome)
    for i in range(2) :
        current_genome, current_barrieres, current_size = mutations(current_genome,current_barrieres, current_size)
        transcriptome = np.array(sim.start_transcribing(INI_file,output_dir))
        current_fitness = fitness(transcriptome,target_profile)
        print(current_fitness, transcriptome)
        print(current_genome)
        hist_fitness.append(current_fitness)
    return []

metropolis()