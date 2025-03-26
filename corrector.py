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
guide_fw = all[14]
guide_rv = all [16]
fw_seq = all[2]
rv_seq = all[4]

A_fw = []; A_rv = []
C_fw = []; C_rv = []
G_fw = []; G_rv = []
T_fw = []; T_rv = []

#### PERHAPS A CONFIDENCE METRIC CAN BE CALCULATED AS THE VALUE OF EACH CHANNEL DIVIDED BY THE TOTAL AND CHOOSE ONWARDS FROM A THRESHOLD OF HIGH REPETITIONS
def confidence(peaks, dol): 
    conl = []
    rm_base = []
    intens = []
    for i in peaks:
        sums = sum([dol[c][i] for c in dol.keys()])
        intens.append(max([dol[c][i] for c in dol.keys()]))
        try:
            conme = max([dol[c][i]/sums for c in dol.keys()])
        except ZeroDivisionError:
            rm_base.append(peaks.index(i))
        conl.append(conme)
    return conl, intens, rm_base

def peak_discovery(dol, guide):
    tmp_lst = list(dol.keys())
    for i in range(0, len(dol.keys())):
        dol[guide[i]] = dol.pop(str(tmp_lst[i]))
    keys = list(dol.keys())
    peaks_key = {}
    lop = []
    for key in keys:
        key_list=[]
        for i in range(2, len(dol[key])-2):
            if dol[key][i] - dol[key][i-1] < 0 and dol[key][i-1] - dol[key][i-2] > 0: #and dol[key][i] != 0 and dol[key][i-1] != 0 and dol[key][i-2] != 0
                key_list.append(i)
        peaks_key[key] = key_list
        lop = lop + peaks_key[key]
    lop.sort()
    return peaks_key, lop

all_peaks = peak_discovery(channels_fw, guide_fw)[1]
new_peaks = all_peaks.copy()
tmp_lst = list(ploc_fw).copy()
for i in all_peaks:
    for j in ploc_fw:
        if i in range(j-5, j+5) and i in new_peaks and j in tmp_lst:
            # print(i, j)
            try:
                tmp_lst.remove(j)
                new_peaks.remove(i)
            except ValueError:
                pass

print(confidence(tmp_lst, channels_fw)[0], confidence(new_peaks, channels_fw)[0])

print(confidence(tmp_lst, channels_fw)[1], confidence(new_peaks, channels_fw)[1])

# for i in range(locs[0][0], locs[0][1]):


# for j in range(locs[1][0], locs[1][1]):


# plt.plot(list(channels_fw.values())[0], "blue")
# plt.show()