######## IMPORTS NEEDED:

import MDAnalysis
from MDAnalysis.analysis.rms import *  #this pulls in an rmsd function
import numpy as np
from MDAnalysis.analysis.align import *
import matplotlib.pyplot as plt
from MDAnalysis.analysis.distances import *

######## GLOBAL SETTINGS:

universe = MDAnalysis.Universe("PATH.prmtop", "PATH.dcd")
np.set_printoptions(threshold=np.inf, precision=2)

prot = universe.select_atoms('protein')

####### PART 1:

#In this part, it is just asked to print the RMSD against the 1st frame

universe.trajectory[0] #first frame
startref = prot.positions
alpha_carbon_1 = prot.select_atoms("name CA").positions

list_carbons = []

for ts in universe.trajectory:
    list_carbons.append(rmsd(alpha_carbon_1, prot.select_atoms("name CA").positions))

np_array_list_carbons = np.array(list_carbons)
#print(np_array_list_carbons)     ##### THIS WILL SHOW YOUR FIRST OUTPUT

####### PART 2:

#In this part, it is asked to show all the RMSDs of All frames

"""universe.trajectory[0] #first frame
alpha_carbon_1 = prot.select_atoms("name CA").positions
universe.trajectory[-1] #last frame
alpha_carbon_2 = prot.select_atoms("name CA").positions

startrmsd = []
endrmsd = []

for ts in universe.trajectory:
    startrmsd.append(rmsd(alpha_carbon_1,prot.select_atoms("name CA").positions))
    endrmsd.append(rmsd(alpha_carbon_2,prot.select_atoms("name CA").positions))

startrmsd_array = np.array(startrmsd)
endrmsd_array = np.array(endrmsd)

print(startrmsd_array)
print(endrmsd_array)"""       ######### THIS IS JUST TO PRINT THE 1st FRAME AND THE LAST FRAME ##########

universe.trajectory[0] #first frame

array_alphas = []

for ts in range(len(universe.trajectory)):
    universe.trajectory[ts]
    alpha_carbon = prot.select_atoms("name CA").positions
    for ts in universe.trajectory:
        array_alphas.append(rmsd(alpha_carbon, prot.select_atoms("name CA").positions))
        np_array_list = np.array(array_alphas)
    #print(np_array_list)           ##### THIS WILL SHOW YOUR SECOND OUTPUT
    array_alphas = []


####### PART 3:

# In this part we give a cutoff and it shows the frame with the highest number of RMSDs within the output (we are going to prepare this part for the PART 4 because it is really similar)

cutoff = 3.3   # You can make this a bit more fancy by creating an input cuttoff    "cutoff = float(input("Write a cutoff number: "))

np_array_list_of_lists = []

for ts in range(len(universe.trajectory)):
    universe.trajectory[ts]
    alpha_carbon = prot.select_atoms("name CA").positions
    for ts in universe.trajectory:
        array_alphas.append(rmsd(alpha_carbon, prot.select_atoms("name CA").positions))
        np_array_list = np.array(array_alphas)
    np_array_list_of_lists.append(np_array_list)                   ##### This is a bit nasty but you can use it to get a really handy list of lists with all the frames
    array_alphas = []

plt.matshow(np_array_list_of_lists)           ##### THIS IS TO PRINT THE WHOLE MATRIX AS A HEATMAP
plt.colorbar()
plt.xlabel("1ENH")
plt.show()

counter_frame = 0
dictionary_frames = {}

for num in range(len(universe.trajectory)):
    dictionary_frames[num] = 0

for list in np_array_list_of_lists:
    for value in list:
        if value<=cutoff:
            dictionary_frames[counter_frame]+=1
        counter_frame+=1
    counter_frame = 0

max_value = 0
frame = 0

for key in dictionary_frames.keys():
    if dictionary_frames[key] > max_value:
        max_value = dictionary_frames[key]
        frame = key

#print("Frame " + str(frame) + " is within " + str(cutoff) + " of " + str(max_value) + " frames")    ##### THIS WILL SHOW YOUR THIRD OUTPUT

####### PART 4:

# Similar to Part 3, but making a small list

cutoff = float(3.3)   # You can make this a bit more fancy by creating an input cuttoff    "cutoff = float(input("Write a cutoff number: "))

np_array_list_of_lists = []

for ts in range(len(universe.trajectory)):
    universe.trajectory[ts]
    alpha_carbon = prot.select_atoms("name CA").positions
    for ts in universe.trajectory:
        array_alphas.append(rmsd(alpha_carbon, prot.select_atoms("name CA").positions))
        np_array_list = np.array(array_alphas)
    np_array_list_of_lists.append(np_array_list)                   ##### This is a bit nasty but you can use it to get a really handy list of lists with all the frames
    array_alphas = []

counter_frame = 0
dictionary_frames = {}

for num in range(len(universe.trajectory)):
    dictionary_frames[num] = 0

for list in np_array_list_of_lists:
    for value in list:
        if value<=cutoff:
            dictionary_frames[counter_frame]+=1
        counter_frame+=1
    counter_frame = 0

max_value = 0
frames = 0
list_of_higher_values = {}

for n in range(10):
    for key in dictionary_frames.keys():
        if dictionary_frames[key] > max_value and key not in list_of_higher_values.keys():
            max_value = dictionary_frames[key]
            frame = key
    list_of_higher_values[frame] = max_value
    max_value = 0

#for key in list_of_higher_values.keys():
    #print("Frame " + str(key) + " is within " + str(cutoff) + " of " + str(dictionary_frames[key]) + " frames")    ##### THIS WILL SHOW YOUR FOURTH OUTPUT

"""for key in dictionary_frames.keys():
    print("Frame " + str(key) + " is within " + str(cutoff) + " of " + str(dictionary_frames[key]) + " frames")"""   #### THIS WOULD BE TO PRINT THE WHOLE LIST

####### PART 5:

# This part is takes a the previous parts (part 4) but it needs to implement an algorithm to get the different features only, considering the range of the cutoff:

cutoff = 1.5   # You can make this a bit more fancy by creating an input cuttoff    "cutoff = float(input("Write a cutoff number: "))

np_array_list_of_lists = []

for ts in range(len(universe.trajectory)):
    universe.trajectory[ts]
    alpha_carbon = prot.select_atoms("name CA").positions
    for ts in universe.trajectory:
        array_alphas.append(rmsd(alpha_carbon, prot.select_atoms("name CA").positions))
        np_array_list = np.array(array_alphas)
    np_array_list_of_lists.append(np_array_list)                   ##### This is a bit nasty but you can use it to get a really handy list of lists with all the frames
    array_alphas = []

counter_frame = 0
counter_frame2 = 0
dictionary_frames = {}

for num in range(len(universe.trajectory)):
    dictionary_frames[num] = []

for list in np_array_list_of_lists:
    for value in list:
        if value<=cutoff:
            dictionary_frames[counter_frame].append(counter_frame2)
        counter_frame+=1
    counter_frame = 0
    counter_frame2 += 1

counter_frame2=0

dictionary_frames_non_overlapped = {}

for key in dictionary_frames.keys():
    if dictionary_frames[key][0] not in dictionary_frames_non_overlapped.keys():
        dictionary_frames_non_overlapped[dictionary_frames[key][0]]=1
    elif dictionary_frames[key][0] in dictionary_frames_non_overlapped.keys():
        dictionary_frames_non_overlapped[dictionary_frames[key][0]]+=1

"""for key in dictionary_frames_non_overlapped.keys():
    print("Frame " + str(key) + " is within " + str(cutoff) + " of " + str(dictionary_frames_non_overlapped[key]) + " frames")"""



