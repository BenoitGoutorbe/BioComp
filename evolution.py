import simulation as sim
import os
import numpy as np
import random
import time

INI_file = "params.ini" #sys.argv[1]
output_dir = "output" #sys.argv[2]
trial_dir = "tousgenesidentiques" #sys.argv[2]

config = sim.read_config_file(INI_file)
TSS_file = config.get('INPUTS','TSS')
TTS_file = config.get('INPUTS','TTS')
GFF_file = config.get('INPUTS', 'GFF')
Prot_file = config.get('INPUTS', 'BARR_FIX')

TSS_file_init = config.get('INPUTS','TSS_init')
TTS_file_init = config.get('INPUTS','TTS_init')
GFF_file_init = config.get('INPUTS', 'GFF_init')
Prot_file_init = config.get('INPUTS', 'BARR_FIX_init')

P_DEL = 0.0
P_INS = 0.0
P_INV = 1.0
SIZE_INDEL = 1000

TOTAL_TIME = time.time()
COUNT_TIME = time.time()

# reinitialiser tousgenesidentiques
def init_tousgenesidentiques():
     #lire les fichiers de tousgenesidentiques_init
    print('REINITIALISATION DU GENOME')
    tss = sim.load_gff(TSS_file_init)
    tts = sim.load_tab_file(TTS_file_init)
    pos_start = tss['TSS_pos'].values
    pos_end = tts['TTS_pos'].values
    orientation = sim.str2num(tss['TUorient'].values)
    positions_genes = [[pos_start[i],pos_end[i],orientation[i]] for i in range(len(pos_start))]
    gff_df_raw = sim.load_gff(GFF_file_init)
    size_genome = int(list(gff_df_raw)[4])
    prot = sim.load_tab_file(Prot_file_init)
    positions_barrieres = prot['prot_pos'].values
    
    #écrire les fichiers de tousgenesidentiques (TSS, TTS et GFF)   
    write_positions(positions_genes, positions_barrieres, size_genome)
    return []

# retourne la fitness du profil en le comparant au profil cible
def fitness (profile, target_profile):
    #difference entre profil et profil cible
    dif = profile/np.sum(profile) - target_profile
    # fitness = 1/distance
    return 1/np.sum(np.multiply(dif,dif))


#lire les fichiers de tousgenesidentiques
def read_positions () :
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

#écrire tous les fichiers du dossier tousgenesidentiques (TSS, TTS et GFF)
def write_positions (pos_genes, pos_bar, size_genome) :
    write_TTS(pos_genes)
    write_TSS(pos_genes)
    write_tousgenesidentiques(pos_genes, size_genome)
    write_prot(pos_bar)
    return []

# ecrire les genes dans TTS
def write_TTS(positions_genes) :
    tts = sim.load_tab_file(TTS_file)
    prob = tts['TTS_proba_off'].values
    file = open(os.path.join(trial_dir,'TTS.dat'), 'w+')
    file.write("TUindex\tTUorient\tTTS_pos\tTTS_proba_off\n")
    for i, v in enumerate(positions_genes):
        file.write(str(i)+ "\t")
        if(v[2] == +1):
            file.write("+\t")
        else:
            file.write("-\t")
        file.write(str(v[1]) + "\t" + str(prob[i]) + "\n")
    file.close()
    return []

# ecrire les genes dans TSS
def write_TSS(positions_genes) :
    tss = sim.load_gff(TSS_file)
    strength = tss['TSS_strength'].values
    file = open(os.path.join(trial_dir,'TSS.dat'), 'w+')
    file.write("TUindex\tTUorient\tTSS_pos\tTSS_strength\n")
    for i, v in enumerate(positions_genes):
        file.write(str(i)+ "\t")
        if(v[2] == +1):
            file.write("+\t")
        else:
            file.write("-\t")
        file.write(str(v[0]) + "\t" + str(strength[i]) + "\n")
    file.close()
    return []

# ecrire les genes dans tousgenesidentiques.gff
def write_tousgenesidentiques(positions_genes, size_genome) :
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

# ecrire les genes dans prot
def write_prot(posbar) :
    file = open(os.path.join(trial_dir,'prot.dat'), 'w+')
    file.write("prot_name\tprot_pos\n")
    for i,v in enumerate(posbar):
        file.write("hns\t" + str(v) + "\n")
    file.close()
    return []

# effectue une mutation sur le genome
# inversion avec une probabilite P_INV
# deletion avec une probabilite P_DEL
# insertion avec une probabilite P_INS
def mutations(pos_genes, pos_barrieres, size_genome) :
    reset_ctime()
    
    positions_genes= list(pos_genes)
    positions_barrieres = list(pos_barrieres)
    # detection des sites codants et non codants
    codant = []
    for g in positions_genes :
        if g[2] > 0 : # si le gene est dans le sens positif
            codant = codant + list(range(g[0],g[1]+1))
        else : # si le gene est dans le sens negatif
            codant = codant + list(range(g[0],g[1]-1, -1))
    non_codant = [i for i in range(1,size_genome+1) if i not in codant]
    
    #inversion
    if random.random() < P_INV :
        # choix de deux bornes dans l'ordre croissant dans le non codant
        bornes = random.sample(non_codant,2)
        bornes.sort()
        print('\nXX INVERSION SUR L\'INTERVALLE', bornes, " XX")
        # inversion des genes
        for g in positions_genes :
            if g[0] > bornes[0] and g[1] < bornes[1]:
                g[0] = bornes[1] - (g[0] - bornes[0])
                g[1] = bornes[1] - (g[1] - bornes[0])
                g[2] = -1.0 * g[2]
        # inversion des barrieres
        for i in range(len(positions_barrieres)):
            if positions_barrieres[i] > bornes[0] and positions_barrieres[i] < bornes[1] :
                positions_barrieres[i] = bornes[1] - (positions_barrieres[i] - bornes[0])
    
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
        for i in range(len(positions_barrieres)):
            if positions_barrieres[i] > position_del :
                positions_barrieres[i] = positions_barrieres[i] - SIZE_INDEL
    #insertion
    if random.random() < P_INS:
        position_ins = random.choice(non_codant)
        print('INSERTION', position_ins)
        size_genome = size_genome + SIZE_INDEL
        for g in positions_genes:
            if g[0] > position_ins and g[1] > position_ins:
                g[0] = g[0] + SIZE_INDEL
                g[1] = g[1] + SIZE_INDEL
        for i in range(len(positions_barrieres)):
            if positions_barrieres[i] > position_ins:
                positions_barrieres[i] = positions_barrieres[i] + SIZE_INDEL
    write_positions(positions_genes, positions_barrieres, size_genome)
    return positions_genes, positions_barrieres, size_genome
    


# algorithme Metropolis sur le genome
def metropolis() :
    
    reset_ctime()
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
    affichage_genome()
    affichage_expr(transcriptome)
    print('\nFITNESS :', current_fitness)
    for i in range(2) :
        current_genome, current_barrieres, current_size = mutations(current_genome,current_barrieres, current_size)
        transcriptome = np.array(sim.start_transcribing(INI_file,output_dir))
        current_fitness = fitness(transcriptome,target_profile)
        affichage_genome()
        affichage_expr(transcriptome)
        print('\nFITNESS :', current_fitness)
        hist_fitness.append(current_fitness)
    return []

# reinitialise le compteur
def reset_ctime():
    COUNT_TIME = time.time()
    return()

# print le temps du compteur s'il n'y a pas d'argument
# renvoie le temps du compteur sinon 
def ctime(ARG = 0):
    if ARG == 0: 
        print('.. ',round(time.time() - COUNT_TIME,3),'s ..\n')
        return()
    else:
        return(round(time.time() - COUNT_TIME,3))


# print le temps total du programme s'il n'y a pas d'argument
# renvoie le temps total du compteur sinon 
def total_time(ARG = 0):
    if ARG == 0: 
        print('\n-- ',round(time.time() - TOTAL_TIME,1),' --')
        return()
    else:
        return(round(time.time() - TOTAL_TIME,1))

# afficher le genome de maniere schematique avec les genes et les proteines
def affichage_genome():
    GENES = sim.load_gff(GFF_file)
    PROTS = sim.load_tab_file(Prot_file)
    
    barrieres= list(PROTS['prot_pos'])
    
    names = list(GENES.loc[:,'ID=id0;Name=tousgenesidentiques'])
    # reduction des noms pour avoir juste g1, g2, g3, etc.
    for i in range(len(names)):
        a = names[i]
        names[i] = a[11:]
    sites1 = list(GENES.loc[:,'1'])
    sites2 = list(GENES.loc[:,'30000'])
    indices = list(GENES.loc[:,'+'])

    # vecteur de toutes les positions, genes comme proteines
    allPos = sites1 + barrieres
    # distinction entre les genes et les proteines
    length_sites1 = len(sites1)
    
    print('\nGENOM:  ', end = '')
    for i in np.argsort(allPos):
        # si c'est un gene
        if i < length_sites1:
            if indices[i] == "+":
                print('   |',names[i],'>    ', end = '', sep = '')
            else:
                print('   <',names[i],'|    ', end = '', sep = '')
        # si c'est une proteine
        else:
            print('.  ', end = '',sep = '')
            
    print('\nGENES:  ', end = '')
    # print des valeurs des sites
    for i in np.argsort(allPos):
        # si c'est un gene
        if i < length_sites1:
            if indices[i] == "+":
                print("|",sites1[i],'-', sites2[i], '>   ', end = '', sep = '')
            else:
                print("<",sites2[i],'-', sites1[i], '|   ', end = '', sep = '')
                
    print('\nPROTE:  ', end = '')
    # print des valeurs des proteines
    for i in np.argsort(allPos):
        # si c'est une proteine
        if i >= length_sites1:
            print('.',barrieres[i - length_sites1],'.         ', end = '',sep = '')
    print('\n')
    return()

def affichage_expr(transcriptome):
    print('EXPRESSION : [g1  g2  g3  g4  g5  g6  g7  g8  g9  g10]')
    print('            ',transcriptome)
    return()

init_tousgenesidentiques()
metropolis()

