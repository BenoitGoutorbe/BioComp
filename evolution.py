'''
Version 1.0

python3 evolution.py n [nombre d iterations] m [nombre de mesures par iteration]

Le temps sera d'a peu pres n*m*3 secondes
'''

import simulation as sim
import os
import numpy as np
import random
import time
import random
import matplotlib.pyplot as plt
import sys


##################
### PARAMETRES ###
##################


INI_file = "params.ini" #sys.argv[1]
output_dir = "output" #sys.argv[2]

# hyp_dir est le dossier qui stocke le prochain pas hypothetique (metropolis)
# c'est le dossier que lit le programme pour faire la transcrition, cf le code de Sam
hyp_dir = "tousgenesidentiques"
# old_dir est le dossier qui stocke la dernière evolution effective de monsieur tousgenesidentiques
# utile pour avoir deux versions du genomes et pouvoir revenir un pas en arriere
old_dir = "tousgenesidentiques_old" 

config = sim.read_config_file(INI_file)
TSS_file = config.get('INPUTS','TSS')
TTS_file = config.get('INPUTS','TTS')
GFF_file = config.get('INPUTS', 'GFF')
Prot_file = config.get('INPUTS', 'BARR_FIX')

TSS_file_init = config.get('INPUTS','TSS_init')
TTS_file_init = config.get('INPUTS','TTS_init')
GFF_file_init = config.get('INPUTS', 'GFF_init')
Prot_file_init = config.get('INPUTS', 'BARR_FIX_init')

TSS_file_old = config.get('INPUTS','TSS_old')
TTS_file_old = config.get('INPUTS','TTS_old')
GFF_file_old = config.get('INPUTS', 'GFF_old')
Prot_file_old = config.get('INPUTS', 'BARR_FIX_old')

P_DEL = 0.0
P_INS = 0.0
P_INV = 1.0
SIZE_INDEL = 1000

TOTAL_TIME = time.time() # on l'affiche avec total_time()
COUNT_TIME = time.time() # compteur. On le reset avec reset_ctime() et l'affiche avec ctime()


#####################################
### FONCTIONS ECRITURE ET LECTURE ###
#####################################


# reinitialiser un dossier (hyp_dir ou old_dir) a partir de tousgenesidentiques_init 
def init_tousgenesidentiques(writing_dir):
    print('REINITIALISATION DU DOSSIER', writing_dir)
    
    #lire les fichiers de tousgenesidentiques_init
    tss = sim.load_gff(TSS_file_init)
    tts = sim.load_tab_file(TTS_file_init)
    gff_df_raw = sim.load_gff(GFF_file_init)
    prot = sim.load_tab_file(Prot_file_init)
    
    pos_start = tss['TSS_pos'].values
    pos_end = tts['TTS_pos'].values
    orientation = sim.str2num(tss['TUorient'].values)
    
    positions_genes = [[pos_start[i],pos_end[i],orientation[i]] for i in range(len(pos_start))]
    size_genome = int(list(gff_df_raw)[4])
    positions_barrieres = prot['prot_pos'].values
    
    #écrire les fichiers (TSS, TTS et GFF)   
    write_positions(positions_genes, positions_barrieres, size_genome,writing_dir)
    return []

# lire les fichiers de tousgenesidentiques
# retourne positions_genes, positions_barrieres, size_genome
def read_positions (reading_dir) :
    if reading_dir == old_dir:
        tss = sim.load_gff(TSS_file_old)
        tts = sim.load_tab_file(TTS_file_old)
        gff_df_raw = sim.load_gff(GFF_file_old)
        prot = sim.load_tab_file(Prot_file_old)
    else:
        tss = sim.load_gff(TSS_file)
        tts = sim.load_tab_file(TTS_file)
        gff_df_raw = sim.load_gff(GFF_file)
        prot = sim.load_tab_file(Prot_file)
    
    pos_start = tss['TSS_pos'].values
    pos_end = tts['TTS_pos'].values
    orientation = sim.str2num(tss['TUorient'].values)
    
    positions_genes = [[pos_start[i],pos_end[i],orientation[i]] for i in range(len(pos_start))]
    size_genome = int(list(gff_df_raw)[4])
    positions_barrieres = prot['prot_pos'].values
    return (positions_genes, positions_barrieres, size_genome)

# écrire tous les fichiers (TSS, TTS et GFF) dans un dossier au choix (hyp_dir ou old_dir)
def write_positions (pos_genes, pos_bar, size_genome, writing_dir) :
    write_TTS(pos_genes,writing_dir)
    write_TSS(pos_genes,writing_dir)
    write_tousgenesidentiques(pos_genes, size_genome,writing_dir)
    write_prot(pos_bar,writing_dir)
    return []

# ecrire les genes dans TTS
def write_TTS(positions_genes, writing_dir) :
    file = open(os.path.join(writing_dir,'TTS.dat'), 'w+')
    file.write("TUindex\tTUorient\tTTS_pos\tTTS_proba_off\n")
    for i, v in enumerate(positions_genes):
        file.write(str(i)+ "\t")
        if(v[2] == +1):
            file.write("+\t")
        else:
            file.write("-\t")
        file.write(str(v[1]) + "\t" + "1.\n")
    file.close()
    return []

# ecrire les genes dans TSS
def write_TSS(positions_genes,writing_dir) :
    tss = sim.load_gff(TSS_file)
    strength = tss['TSS_strength'].values
    file = open(os.path.join(writing_dir,'TSS.dat'), 'w+')
    file.write("TUindex\tTUorient\tTSS_pos\tTSS_strength\n")
    for i, v in enumerate(positions_genes):
        file.write(str(i)+ "\t")
        if(v[2] == +1):
            file.write("+\t")
        else:
            file.write("-\t")
        file.write(str(v[0]) + "\t" + str(strength[i])[1:] + "\n")
    file.close()
    return []

# ecrire les genes dans tousgenesidentiques.gff
def write_tousgenesidentiques(positions_genes, size_genome,writing_dir) :
    file = open(os.path.join(writing_dir,'tousgenesidentiques.gff'), 'w+')
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
def write_prot(posbar, writing_dir) :
    file = open(os.path.join(writing_dir,'prot.dat'), 'w+')
    file.write("prot_name\tprot_pos\n")
    for i,v in enumerate(posbar):
        file.write("hns\t" + str(v) + "\n")
    file.close()
    return []


##################################
### FONCTIONS POUR L EVOLUTION ###
##################################


# retourne la fitness du profil en le comparant au profil cible
def fitness (profile, target_profile):
    # difference entre profil et profil cible
    dif = profile/np.sum(profile) - target_profile
    # fitness = 1/distance
    return 1/np.sum(np.multiply(dif,dif)) # ou 1-np.sum(np.multiply(dif,dif))


# effectue une mesure de fitness m fois
# renvoie un tableau de taille m    
def fitness_sample(genome, barrieres, size, target, m):
    V = []
    write_positions(genome, barrieres, size, hyp_dir)
    for i in range(m):
        transcriptome = np.array(sim.start_transcribing(INI_file, output_dir))
        V.append(fitness(transcriptome, target))
    return(V)


# effectue une mutation sur le genome
# inversion, deletion ou insertion avec des probabilites respectives de P_INV, P_DEL et P_INS
def mutations(pos_genes, pos_barrieres, size_genome) :
    positions_genes = list(pos_genes)
    positions_barrieres = list(pos_barrieres)
    # detection des sites codants
    codant = []
    for g in positions_genes :
        if g[2] > 0 : # si le gene est dans le sens positif
            codant = codant + list(range(g[0],g[1]+1))
        else : # si le gene est dans le sens negatif
            codant = codant + list(range(g[0],g[1]-1, -1))
    
    # detection des sites non-codants
        # ancienne version : non_codant = [i for i in range(1,size_genome+1) if i not in codant]
    sites = range(1,size_genome+1)
    non_codant = list(set(sites)- set(codant))    
    
    #inversion
    if random.random() < P_INV :
        # choix de deux bornes dans l'ordre croissant dans le non codant
        bornes = random.sample(non_codant,2)
        bornes.sort()
        print('\nX INVERSION SUR L\'INTERVALLE', bornes, "X\n")
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
    return positions_genes, positions_barrieres, size_genome
    
# importe le profil cible du fichier environnement
def target():
    envir_file = open(os.path.join('tousgenesidentiques', 'environment.dat'))
    target_profile = []
    for line in envir_file.readlines():
        target_profile.append(float(line.split()[1]))
    target_profile = np.array(target_profile)
    envir_file.close()
    return(target_profile)

# algorithme Metropolis sur le genome
# n est le nombre d'iteration
# m est le nombre de mesures de fitness par iteration
def metropolis(n,m) :
    
    print('\nDebut de simulation. ', n, ' iterations et ', m, 'mesures de fitness par iteration.')
    print('Temps estime :', n*m*3, 'secondes.\n')
    
    # importation du profil cible
    target_profile = target()
    
    # lecture du genome initial
    ini_genome, ini_barrieres, ini_size = read_positions(old_dir)
    
    # mesure de la fitness initiale
    hyp_Vfitness = fitness_sample(ini_genome, ini_barrieres, ini_size, target_profile, m)
    
        # ancienne version (une seule mesure par iteration)
        #~ write_positions(ini_genome, ini_barrieres, ini_size, hyp_dir)
        #~ ini_transcriptome = np.array(sim.start_transcribing(INI_file, output_dir))
        #~ ini_fitness = fitness(ini_transcriptome, target_profile)
    
    # hist_fitness contiendra les fitness *moyennes* au fur et a mesure de la simulation.
    # fitness tab contient toutes les m fitness des n iterations. Soit un tableau n * m
    hist_fitness = [np.mean(hyp_Vfitness)]
    fitnessTab = [hyp_Vfitness]
    
    affichage_genome(old_dir)
    affichage_fitness(hist_fitness[-1])
    
    # ces variables correspondent au genome qui a eu la meilleure fitness moyenne de la simulation
    best_genome = ini_genome
    best_barrieres = ini_barrieres
    best_size = ini_size
    best_fitness = hist_fitness[-1]
    
    for i in range(n) :
        print('\n_ _ _ _ _ Iteration :', i+1, 'sur', n)
        print('_ _ _ _ _', ctime(1), 's')
        
        # recuperation des dernieres valeurs du genome
        old_genome, old_barrieres, old_size = read_positions(old_dir)
        # la mutation de ce genome renvoie le prochain genome hypothetique
        hyp_genome, hyp_barrieres, hyp_size = mutations(old_genome, old_barrieres, old_size)
        
            # ancienne version (une seule mesure par iteration)
            #~ write_positions(hyp_genome, hyp_barrieres, hyp_size, hyp_dir)
            #~ hyp_transcriptome = np.array(sim.start_transcribing(INI_file,output_dir))
            #~ hyp_fitness = fitness(hyp_transcriptome,target_profile)
        
        # Mesures de m valeurs de fitness
        hyp_Vfitness = fitness_sample(hyp_genome, hyp_barrieres, hyp_size, target_profile, m)
        # Fitness moyenne sur les m mesures
        mean_hyp_Vfitness = np.mean(hyp_Vfitness)
        
        # Comparaison de la fitness moyenne a la fitness precedente
        alpha = mean_hyp_Vfitness / hist_fitness[-1]
        
        affichage_genome(hyp_dir)
        affichage_fitness(mean_hyp_Vfitness)
        
        # alpha > 1 : la fitness a augmente
        # il y a aussi une chance alpha/10 de garder le genome meme s'il est moins bon
        if alpha > 1 or random.random() < alpha/10 :
            
            fitnessTab.append(hyp_Vfitness)
            hist_fitness.append(mean_hyp_Vfitness)
            
            # le genome hypothetique devient effectif : on l'ecrit dans old_dir
            write_positions(hyp_genome, hyp_barrieres, hyp_size, old_dir)
            
            # si c'est la meilleure fitness moyenne, on garde en stock le genome
            if mean_hyp_Vfitness > best_fitness :
                best_genome = hyp_genome
                best_barrieres = hyp_barrieres
                best_size = hyp_size
                best_fitness = mean_hyp_Vfitness
            
            if alpha > 1:
                print('\nAugmentation de la fitness moyenne, mutation retenue')
            else:
                print('\nDiminution de la fitness moyenne, mais mutation retenue')
        else:
            print('\nDiminution de la fitness moyenne, mutation non retenue')
    
    
    print('\n\n_ _ _ _ _ Le meilleur genome')
    
    # ecriture du dernier genome dans le dossier hyp pour pouvoir l'afficher
    write_positions(best_genome, best_barrieres, best_size, hyp_dir)
    affichage_genome(hyp_dir)
    affichage_fitness(best_fitness)
    
    
    print('\nFin de la simulation.\nTemps total :',ctime(1),'s')
    
    # plot de la simulation complete
    affichage_simulation(fitnessTab,hist_fitness)
    
    return []


###############################
### FONCTIONS POUR LE TEMPS ###
###############################


# reinitialise le compteur
def reset_ctime():
    COUNT_TIME = time.time()
    return()

# print le temps du compteur s'il n'y a pas d'argument
# renvoie le temps du compteur sinon 
def ctime(ARG = 0):
    if ARG == 0:
        print('..',round(time.time() - COUNT_TIME,3),'s ..')
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


##################################
### FONCTIONS POUR L AFFICHAGE ###
##################################


# afficher le genome de maniere schematique avec les genes et les proteines
def affichage_genome(reading_dir):
    if reading_dir == old_dir:
        GENES = sim.load_gff(GFF_file_old)
        PROTS = sim.load_tab_file(Prot_file_old)
    else:
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
    print('')
    return()


# afficher l'expression
def affichage_expr(transcriptome):
    print('EXPRESSION : [g1  g2  g3  g4  g5  g6  g7  g8  g9  g10]')
    print('            ',transcriptome)
    return()

# afficher la fitness
def affichage_fitness(fitness):
    print('\nFITNESS :', fitness)
    return()


# affichage du total de la simulation
def affichage_simulation(fitnessTab,hist_fitness):
    axes = plt.gca()
    xL = len(hist_fitness)-1    
    axes.set_xlim([-0.05*xL,xL*1.05])
    axes.set_ylim([0,np.max(fitnessTab)*1.05])
    plt.ylabel('Mesures de fitness et moyenne')
    plt.xlabel('Iteration')
    
    plt.plot(hist_fitness)
    for i,Vfit in enumerate(fitnessTab) :
        plt.scatter(np.ones(len(Vfit))*(i),Vfit)
    
    print('Plot ouvert dans une nouvelle fenetre. L\'enregistrer avant de fermer la fenetre.')
    
    plt.show()
    return()


######################
### FONCTIONS MAIN ###
######################


init_tousgenesidentiques(hyp_dir)
init_tousgenesidentiques(old_dir)
metropolis(int(sys.argv[1]),int(sys.argv[2]))


