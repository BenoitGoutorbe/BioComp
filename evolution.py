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

P_INV = 0.5
P_DEL = 0.25
P_INS = 0.25

SIZE_INDEL = 500

TOTAL_TIME = time.time() # on l'affiche avec total_time()
COUNT_TIME = time.time() # compteur. On le reset avec reset_ctime() et l'affiche avec ctime()

SIMUL_NAME = time.strftime("J%d_%Hh_%Mmn")

# ligne pour la sauvegarde des infos dans des fichiers .txt
result_dir = "output_" + SIMUL_NAME
os.mkdir(result_dir)

file_fit = open(os.path.join(result_dir,'HIST_fitness.txt'), 'w+')
file_console = open(os.path.join(result_dir,'HIST_console.txt'), 'w+')

#####################################
### FONCTIONS ECRITURE ET LECTURE ###
#####################################


# reinitialiser un dossier (hyp_dir ou old_dir) a partir de tousgenesidentiques_init 
def init_tousgenesidentiques(writing_dir):
    print_save(str('REINITIALISATION DU DOSSIER '+ writing_dir), file_console )
    
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
    #return np.exp(-np.sum(np.multiply(dif,dif)))
    pente = 10 # une pente de 10 permet d'avoir une fitness qui se repartit generalement entre 0 et 1
    return 2 / (1 + np.exp( pente * np.sqrt( np.sum( np.multiply(dif,dif) ) ) ) )

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
    
    # choix aleatoire du type de mutation selon P_INV, P_DEL, P_INS
    p = random.random()
    # inversion
    if 0 <= p and p < P_INV :
        type_mut = "INV"
        # choix de deux bornes dans l'ordre croissant dans le non codant
        bornes = random.sample(non_codant,2)
        bornes.sort()
        print_save(str('\n[x] INVERSION SUR L\'INTERVALLE '+ str(bornes) + " [x]\n"), file_console )
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
    
    # deletion
    if P_INV <= p and p < P_INV + P_DEL:
        type_mut = "DEL"
        # choix d'un site de taille SIZE_INDEL a retirer, qui ne soit pas en zone codante
        # au bout de 1000 essais au compteur, on arrete d essayer
        compteur = 1
        CORRECT_DEL = False
        while CORRECT_DEL == False and compteur <= 1000 :
            CORRECT_DEL = True
            position_del = random.choice(non_codant)
            # prise en compte des genes
            if(position_del + SIZE_INDEL in codant):
                CORRECT_DEL = False
            # prise en compte des barrieres
            for bar in positions_barrieres:
                if position_del <= bar and bar <= position_del + SIZE_INDEL:
                    CORRECT_DEL = False
            
            compteur += 1
        
        if compteur <= 1000:
            print_save('\n[-] DELETION DE ' + str(SIZE_INDEL) + ' SITES DE ' + str(position_del) + ' A ' + str(position_del + SIZE_INDEL) +' [-]\n', file_console )
            size_genome = size_genome - SIZE_INDEL
            for g in positions_genes:
                if g[0] > position_del and g[1] > position_del:
                    g[0] = g[0] - SIZE_INDEL
                    g[1] = g[1] - SIZE_INDEL
            
            for i in range(len(positions_barrieres)):
                if positions_barrieres[i] > position_del :
                    positions_barrieres[i] = positions_barrieres[i] - SIZE_INDEL
        else:
            print_save('\nPas de deletion possible car aucune region de taille ' + str(SIZE_INDEL) + ' n\'a été trouvée.', file_console )
    
    # insertion
    if P_INV + P_DEL <= p and p < P_INV + P_DEL + P_INS:
        type_mut = "INS"
        # choix d'un site pour l insertion
        position_ins = random.choice(non_codant)
        print_save('\n[+] INSERTION DE ' + str(SIZE_INDEL) + ' SITES EN POSITION ' + str(position_ins) + ' [+]\n', file_console )
        size_genome = size_genome + SIZE_INDEL
        for g in positions_genes:
            if g[0] > position_ins and g[1] > position_ins:
                g[0] = g[0] + SIZE_INDEL
                g[1] = g[1] + SIZE_INDEL
        
        
        for i in range(len(positions_barrieres)):
            if positions_barrieres[i] > position_ins:
                positions_barrieres[i] = positions_barrieres[i] + SIZE_INDEL
    return positions_genes, positions_barrieres, size_genome , type_mut
    
# importe le profil cible du fichier environnement
def target():
    envir_file = open(os.path.join('tousgenesidentiques', 'environment.dat'))
    target_profile = []
    for line in envir_file.readlines():
        target_profile.append(float(line.split()[1]))
    target_profile = np.array(target_profile)
    envir_file.close()
    return(target_profile)

def print_save(str1, file1 = 0):
    if file1 != 0:
        file1.write(str1 + "\n")
    print(str1)
    
    return()
# algorithme Metropolis sur le genome
# n est le nombre d'iteration
# m est le nombre de mesures de fitness par iteration
def metropolis(n,m) :
    
    # initialisation des fichiers d ecriture
    file_fit.write("id")
    for i in range(m):
        file_fit.write("\tfit" + str(i))
        
    file_fit.write("\tmeanFit\tret\n")

    
    print_save('\nDebut de simulation. ' + str(n) + ' iterations et ' + str(m) + ' mesures de fitness par iteration.', file_console)
    print_save('Temps estime : ' + str(n*m*3) + ' secondes.\n', file_console)
    
    
    
    
    
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
    
    # Comptage des succes de chaque type de mutation
    MUT_NUM = np.zeros((1,3))
    MUT_SUCCESS = np.zeros((1,3))
    DICO = dict()
    DICO["INV"] = 0
    DICO["DEL"] = 1
    DICO["INS"] = 2
    
    for i in range(n) :
        # print de l'id
        file_fit.write(str(i+1))
        
        print_save('\n_ _ _ _ _ Iteration : ' + str(i+1) + ' sur ' + str(n), file_console)
        print_save('_ _ _ _ _' + str(ctime(1)) + ' s ', file_console)
        
        # recuperation des dernieres valeurs du genome
        old_genome, old_barrieres, old_size = read_positions(old_dir)
        # la mutation de ce genome renvoie le prochain genome hypothetique
        hyp_genome, hyp_barrieres, hyp_size, type_mut = mutations(old_genome, old_barrieres, old_size)
        
            # ancienne version (une seule mesure par iteration)
            #~ write_positions(hyp_genome, hyp_barrieres, hyp_size, hyp_dir)
            #~ hyp_transcriptome = np.array(sim.start_transcribing(INI_file,output_dir))
            #~ hyp_fitness = fitness(hyp_transcriptome,target_profile)
        
        # Mesures de m valeurs de fitness
        hyp_Vfitness = fitness_sample(hyp_genome, hyp_barrieres, hyp_size, target_profile, m)
        
        for fit in hyp_Vfitness:
            file_fit.write("\t" + str(fit))
        
        MUT_NUM[0,DICO[type_mut]] += 1
        # Fitness moyenne sur les m mesures
        mean_hyp_Vfitness = np.mean(hyp_Vfitness)
        
        # ecriture de la fitness moyenne
        file_fit.write("\t" + str(mean_hyp_Vfitness))
        
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
                print_save('\nAugmentation de la fitness moyenne, mutation retenue', file_console)
                MUT_SUCCESS[0,DICO[type_mut]] = 1
                file_fit.write('\t2\n')
            else:
                print_save('\nDiminution de la fitness moyenne, mais mutation retenue', file_console)
                file_fit.write('\t1\n')
        else:
            print_save('\nDiminution de la fitness moyenne, mutation non retenue', file_console)
            file_fit.write('\t0\n')
    
    
    print_save('\n\n_ _ _ _ _ Le meilleur genome', file_console )
    
    # ecriture du best genome dans le dossier hyp afin de l'afficher
    write_positions(best_genome, best_barrieres, best_size, hyp_dir)
    affichage_genome(hyp_dir)
    affichage_fitness(best_fitness)
    
    print_save('\n Nombre de INV, DEL et INS et pourcentages de succes respectifs :', file_console )
    RES = MUT_SUCCESS/MUT_NUM
    print_save(str(MUT_NUM[0,:]), file_console )     
    print_save(str(RES[0,:]), file_console )    
    
    print_save('\nFin de la simulation.\nTemps total : ' + str(ctime(1)) + ' s', file_console )
    
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
        print_save('.. ' + str(round(time.time() - COUNT_TIME,3) ) + ' s ..', file_console )
        return()
    else:
        return(round(time.time() - COUNT_TIME,3))


# print le temps total du programme s'il n'y a pas d'argument
# renvoie le temps total du compteur sinon 
def total_time(ARG = 0):
    if ARG == 0: 
        print_save('\n-- ' + str(round(time.time() - TOTAL_TIME,1) )  + ' --', file_console )
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
    gSIZE = GENES.axes[1][4]
    sites2 = list(GENES.loc[:,gSIZE])
    indices = list(GENES.loc[:,'+'])
    
    # vecteur de toutes les positions, genes comme proteines
    allPos = sites1 + barrieres
    # distinction entre les genes et les proteines
    length_sites1 = len(sites1)
    
    STRING_GENOM = '\nGENOM:  '
    for i in np.argsort(allPos):
        # si c'est un gene
        if i < length_sites1:
            if indices[i] == "+":
                STRING_GENOM += '   |' + str(names[i]) + '>    '
            else:
                STRING_GENOM += '   <' + str(names[i]) + '|    '
        # si c'est une proteine
        else:
            STRING_GENOM += '.  '
    print_save(STRING_GENOM, file_console)
    
    STRING_GENES = 'GENES:  '
    # print des valeurs des sites
    for i in np.argsort(allPos):
        # si c'est un gene
        if i < length_sites1:
            if indices[i] == "+":
                STRING_GENES += "|" + str(sites1[i]) + '-' + str(sites2[i]) + '>   '
            else:
                STRING_GENES += "<" + str(sites2[i]) + '-' + str(sites1[i]) + '|   '
    print_save(STRING_GENES, file_console)
    
    STRING_PROTE = 'PROTE:  '
    # print des valeurs des proteines
    for i in np.argsort(allPos):
        # si c'est une proteine
        if i >= length_sites1:
            STRING_PROTE += '.' + str(barrieres[i - length_sites1]) + '.         '
    print_save(STRING_PROTE, file_console)
    
    print_save('', file_console)
    return()


# afficher l'expression
def affichage_expr(transcriptome):
    print_save('EXPRESSION : [g1  g2  g3  g4  g5  g6  g7  g8  g9  g10]', file_console )
    print_save('             ' + str(transcriptome), file_console )
    return()

# afficher la fitness
def affichage_fitness(fitness):
    print_save('\nFITNESS : ' + str(fitness), file_console )
    return()


# affichage du total de la simulation
def affichage_simulation(fitnessTab,hist_fitness):
    axes = plt.gca()
    xL = len(hist_fitness)-1    
    axes.set_xlim([-0.05*xL,xL*1.05])   
    axes.set_ylim([0,np.max(fitnessTab)*1.05])
    plt.ylabel('Mesures de fitness et moyenne')
    plt.xlabel('Iteration')
    plt.title('Evolution de la fitness')
    plt.plot(hist_fitness)
    
    for i,Vfit in enumerate(fitnessTab) :
        plt.scatter(np.ones(len(Vfit))*(i),Vfit)
    
    print_save('Graphe ouvert dans une nouvelle fenetre et enregistré avec les résultats en ' + str(result_dir) + ".", file_console )
    
    plt.savefig(result_dir + '/GRAPHE_' + time.strftime("J%d_%Hh_%Mmn") + ".png")
    plt.show()
    return()


######################
### FONCTIONS MAIN ###
######################


init_tousgenesidentiques(hyp_dir)
init_tousgenesidentiques(old_dir)
metropolis(int(sys.argv[1]),int(sys.argv[2]))


file_fit.close()
file_console.close()
