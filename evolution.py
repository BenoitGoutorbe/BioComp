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

def mutations () :
    #generer les mutations et réécrire les fichiers de tousgenesidentiques
    return []

def metropolis() :
    envir_file = open(os.path.join('tousgenesidentiques', 'environment.dat'))
    target_profile = []
    for line in envir_file.readlines():
        target_profile.append(float(line.split()[1]))
    target_profile = np.array(target_profile)
    envir_file.close()
    mutations()
    #simulation :
    transcriptome = np.array(sim.start_transcribing(INI_file, output_dir))
    current_fitness = fitness(transcriptome, target_profile)
    print(current_fitness)
    return []

metropolis()