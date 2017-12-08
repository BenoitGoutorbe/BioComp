import simulation as sim
import os
import numpy as np

INI_file = "params.ini" #sys.argv[1]
output_dir = "output" #sys.argv[2]


env_file = open(os.path.join('tousgenesidentiques','environment.dat'))
target_profile = []
for line in env_file.readlines():
    target_profile.append(float(line.split()[1]))
target_profile= np.array(target_profile)

def fitness (profile):
    #difference entre profil et profil cible
    return np.mean(np.multiply(profile/np.sum(profile) -target_profile))

def mutations () :
    #generer les mutations et réécrire les fichiers de tousgenesidentiques
    return []

def metropolis() :
    mutations()
    #sim.start_transcribing(INI_file, output_dir)
    profil_cible = []  # lire envirnoment.dat
    return []
