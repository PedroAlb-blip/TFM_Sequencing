#### This file will take the values stored from parser.py and produce the correction 
#### IDEAS:
#### 1.- Look for local peaks and produce a list containing them (all of them)
#### 2.- Cross check the PLOC lists of both forward and reverse
#### These should provide useful data regarding non-aligned regions
#### For aligned regions:
#### 1.- Eliminate insertions by cross checking the local peak list
#### 2.- Wherever there is a mismatch do the sum of channels and verify the nucleotide in either forward or reverse
#### 3.- Alter the forward and reverse sequence and make a new alignment
import ast
import matplotlib.pyplot as plt
file = open("Sp101b-VP7.txt", "r")
data = file.read()
all = list(filter(('').__ne__, data.split('\n')))
###### IMPORTANT INFORMATION IN: 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22 and 24
channels = ["DATA9", "DATA10", "DATA11", "DATA12"]
channels_fw = ast.literal_eval(all[10])
channels_rv = ast.literal_eval(all[12])
ploc_fw = ast.literal_eval(all[6])
ploc_rv = ast.literal_eval(all[8])
locs = ast.literal_eval(all[24])
guide = all[16]
fw_seq = all[2]
rv_seq = all[4]

A_fw = []; A_rv = []
C_fw = []; C_rv = []
G_fw = []; G_rv = []
T_fw = []; T_rv = []

#### PERHAPS A CONFIDENCE METRIC CAN BE CALCULATED AS THE VALUE OF EACH CHANNEL DIVIDED BY THE TOTAL AND CHOOSE ONWARDS FROM A THRESHOLD OF HIGH REPETITIONS

for i in ploc_fw[0:100]:
    sums = sum([channels_fw[c][i] for c in channels])
    cands = max([channels_fw[c][i]/sums for c in channels if sums !=0])

for i in ploc_rv[-100:-1]:
    sums = sum([channels_rv[c][i] for c in channels])
    cands = max([channels_rv[c][i]/sums for c in channels if sums !=0])
    print(cands)

# for i in range(locs[0][0], locs[0][1]):


# for j in range(locs[1][0], locs[1][1]):


# plt.plot(list(channels_fw.values())[0], "blue")
# plt.show()