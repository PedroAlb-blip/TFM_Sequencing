#### This file will take the values stored from parser.py and produce the correction 
#### IDEAS:
#### 1.- Look for local peaks and produce a list containing them (all of them)
#### 2.- Cross check the PLOC lists of both forward and reverse
#### 3.- This will provide two lists:
####    One of unrealised peaks (present in discovery yet absent from the basecaller)
####    One of false peaks (present in basecaller, absent in discovery)
#### These should provide useful data regarding non-aligned regions
#### For aligned regions:
#### 1.- Eliminate insertions by cross checking the local peak list
#### 2.- Wherever there is a mismatch do the sum of channels and verify the nucleotide in either forward or reverse
#### 3.- Alter the forward and reverse sequence and make a new alignment
import ast
import matplotlib.pyplot as plt
import torch
from torch import nn
from torch.utils.data import DataLoader
from torchvision import datasets, transforms

#### PERHAPS A CONFIDENCE METRIC CAN BE CALCULATED AS THE VALUE OF EACH CHANNEL DIVIDED BY THE TOTAL AND CHOOSE ONWARDS FROM A THRESHOLD OF HIGH REPETITIONS
def confidence(peaks, dol): 
    #### This can be a good parameter to train a nn model with 
    conl = []
    rm_base = []
    intens = []
    peak_1 = peaks[0]
    peak_dist = []
    for i in peaks:
        peak_dist.append(i - peak_1)
        peak_1 = i
        sums = sum([dol[c][i] for c in dol.keys()])
        intens.append(max([dol[c][i] for c in dol.keys()]))
        #### The confidence metric can be calculated as the max over the sum of all intensity values
        try:
            conme = max([dol[c][i]/sums for c in dol.keys()])
        except ZeroDivisionError:
            rm_base.append(peaks.index(i))
        conl.append(conme)
    #### Other values which can be useful for training are the intensities, and distance between peaks
    return conl, intens, peak_dist, rm_base

def rename(dol, guide):
    tmp_lst = list(dol.keys())
    for i in range(0, len(dol.keys())):
        dol[guide[i]] = dol.pop(str(tmp_lst[i]))
    return dol

def peak_discovery(dol):
    #### This function works by returning the indices of possible peaks as per the definition of a local maximum
    keys = list(dol.keys())
    peaks_key = {}
    lop = []
    for key in keys:
        key_list=[]
        for i in range(2, len(dol[key])):
            if dol[key][i] - dol[key][i-1] < 0 and dol[key][i-1] - dol[key][i-2] > 0: #### Wherever the derivatives switch signs there must be a local peak (min or max)
                key_list.append(i-1)
        peaks_key[key] = key_list
        lop = lop + peaks_key[key]
    lop.sort()
    return peaks_key, lop #### This returns the peaks as a dictionary and as a list

def jiggler(lop, dol): #### This function will return a list of jiggled peak locations as severance that they are locally peaks
    new = []
    for i in lop:
        local_max = []
        contester_up = []
        contester_dwn = []
        for c in dol.keys():
            local_max.append(dol[c][i])
            contester_up = contester_up + list(dol[c][i : i + 3])
            contester_dwn = contester_dwn + list(dol[c][i - 3 : i])[::-1]
        if max(contester_up) > max(local_max) or max(contester_dwn) > max(local_max):
            if max(contester_up) >= max(contester_dwn):
                new.append(i + contester_up.index(max(contester_up)) % 3)
            elif max(contester_dwn) > max(contester_up):
                new.append(i - contester_dwn.index(max(contester_dwn)) % 3)
        else:
            new.append(i)
    return new

#### This lines of code cross check the discovered peaks and the ploc from the basecaller 
def cross_check(dol, lop):
    all_peaks = peak_discovery(dol)[1]
    new_peaks = all_peaks.copy()
    ploc_fw = jiggler(lop, dol)
    tmp_lst = list(lop.copy())
    full_on = []
    for i in all_peaks:
        for j in ploc_fw:
            if j in list(range(i-2, i+3)) and i in new_peaks and j in tmp_lst:
                try:
                    tmp_lst.remove(j)
                    new_peaks.remove(i)
                    full_on.append(i)
                except ValueError:
                    pass
    return full_on, new_peaks, tmp_lst


def width(dol, lop):
    amplitude = []
    for i in lop:
        peaks_dist = []
        j = i
        for c in dol.keys():
            if dol[c][i] is max([dol[k][i] for k in dol.keys()]):
                peaks_dist.append(i)
                while dol[c][i] > 0:
                    i = i + 1
                while dol[c][j] > 0:
                    j = j - 1
                amplitude.append(i - j)
                break
    return amplitude

def filterer(dol, lop, alignment, locs):  
    better_peaks = jiggler(lop, dol)
    every = confidence(better_peaks, dol)
    conf = every[0]
    intens = every[1]
    peakdis = every[2]
    amp = width(dol, better_peaks)
    checkers = cross_check(dol, better_peaks)[0]
    duplic = []
    der_1 = []
    der_2 = []
    result = []
    if max([dol[k][0] for k in dol.keys()]) < 100:
        loc = locs[1]
    elif max([dol[k][0] for k in dol.keys()]) > 600:
        loc = locs[0]
    for i in lop:
        if i in checkers:
            duplic.append(0)
        else:
            duplic.append(1)
    for i in better_peaks:
        for c in dol.keys():
            if dol[c][i] == max([dol[k][i] for k in dol.keys()]):
                der_1.append(dol[c][i] - dol[c][i - 1])
                der_2.append(dol[c][i + 1] - dol[c][i])
                break
    j = 0
    for i in range(0,len(lop)):
        if intens[i - j]/conf[i - j] > 1200 or conf[i - j] < 0.6:
            result.append(0)
        elif intens[i - j] < 100:
            result.append(0)
        elif alignment[i + loc[0]] == "*":
            result.append(1)
        else:
            conf.pop(i - j)
            intens.pop(i - j)
            duplic.pop(i - j)
            amp.pop(i - j)
            peakdis.pop(i - j)
            der_1.pop(i - j)
            der_2.pop(i - j)
            j = j + 1
    return conf, intens, duplic, amp, peakdis, der_1, der_2, result
#### True is 0 and False is 1
file = open("Sp101b-VP7.txt", "r")
data = file.read()
all = list(filter(('').__ne__, data.split('\n')))
###### IMPORTANT INFORMATION IN: 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22 and 24
channels_fw = ast.literal_eval(all[10])
channels_rv = ast.literal_eval(all[12])
ploc_fw = ast.literal_eval(all[6])
ploc_rv = ast.literal_eval(all[8])
locs = ast.literal_eval(all[24])
align = all[18]
guide_fw = all[14]
guide_rv = all [16]
fw_seq = all[2]
rv_seq = all[4]


print(filterer(channels_fw, ploc_fw, align, locs)[2])


#### TO TRAIN A NN TO CHECK IF THESE PEAKS EXIST AND THEY CORRELATE TO BASES, THE FOLLOWING PARAEMETERS MUST BE USED:
#### AMPLITUDE OF PEAKS, DISTANCE BETWEEN PEAKS, INTENSITY, CONFIDENCE, DERIVATIVES SIDEWAYS OF PEAK
#### MATCHES WITHIN THE ALIGNMENT CAN BE USED FOR TRANING 

print(locs)


# print(confidence(tmp_lst, channels_fw)[0], confidence(new_peaks, channels_fw)[0])

# print(confidence(tmp_lst, channels_fw)[1], confidence(new_peaks, channels_fw)[1])


# rename(channels_fw, guide_fw)

# print(jiggler(ploc_fw, channels_fw))

# for i in range(locs[0][0], locs[0][1]):


# for j in range(locs[1][0], locs[1][1]):


# plt.plot(list(channels_fw.values())[0], "blue")
# plt.show()