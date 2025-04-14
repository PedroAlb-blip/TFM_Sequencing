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
from tensorflow.python.keras.models import Sequential
from tensorflow.python.keras.layers import Dense, Activation
import numpy as np

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
    newd = {key : [] for key in dol.keys()}
    new = []
    for i in lop:
        for c in dol.keys():
            j = i - 3
            while j < i + 3:
                if dol[c][j] >= dol[c][j+1] and dol[c][j] > 100 and dol[c][j + 1] > 100:
                    newd[c].append(j)
                    new.append(j)
                    break
                else:
                    j = j + 1
    return new

#### This lines of code cross check the discovered peaks and the ploc from the basecaller 
def cross_check(dol, lop):
    all_peaks = peak_discovery(dol)[1]
    new_peaks = all_peaks.copy()
    ploc_fw = jiggler(lop, dol)
    ploc = list(lop).copy()
    full_on = []
    for i in all_peaks:
        for j in ploc_fw:
            if j in list(range(i-2, i+3)) and i in new_peaks and j in ploc:
                try:
                    ploc.remove(j)
                    new_peaks.remove(i)
                    full_on.append(i)
                except ValueError:
                    pass
    return full_on, new_peaks, ploc

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

def filterer(dol, lop, alignment, locs, train = True):  
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
        loc = locs[2] - locs[0][0]
    elif max([dol[k][0] for k in dol.keys()]) > 600:
        loc = locs[2] - locs[1][0]
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
    joined = []
    if train == True:
        for i in range(0,len(lop)):
            joined.append([conf[i - j], intens[i - j], duplic[i - j], amp[i - j], peakdis[i - j], der_1[i - j], der_2[i - j], i - j])
            if intens[i - j] > 1200 or conf[i - j] < 0.6:
                result.append(0)
            elif intens[i - j] < 100:
                result.append(0)
            elif alignment[i + loc] == "*":
                result.append(1)
            else:
                conf.pop(i - j)
                intens.pop(i - j)
                duplic.pop(i - j)
                amp.pop(i - j)
                peakdis.pop(i - j)
                der_1.pop(i - j)
                der_2.pop(i - j)
                joined.pop(-1)
                j = j + 1
        return joined, result
    else:
        for i in range(0,len(lop)):
            joined.append([conf[i], intens[i], duplic[i], amp[i], peakdis[i], der_1[i], der_2[i], i])
        return joined
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


train = filterer(channels_fw, ploc_fw, align, locs)
print(train)
trainX = np.array(train[0])
trainY = np.array(train[1])
model = Sequential()

model.add(Dense(8, input_dim=8, activation='relu'))
model.add(Dense(1, activation='sigmoid'))
model.compile(loss='mean_squared_error', optimizer='adam')
model.fit(trainX, trainY, epochs=200, batch_size=32, verbose=1)

dataPrediction = model.predict(np.array(filterer(channels_fw, ploc_fw, align, locs, train=False)))
pred = list(np.ndarray.flatten(dataPrediction))
dataPrediction = model.predict(np.array(filterer(channels_fw, cross_check(channels_fw, ploc_fw)[1], align, locs, train=False)))
other = list(np.ndarray.flatten(dataPrediction))
print(pred)
print('\n\n')
print(other)
#### TO TRAIN A NN TO CHECK IF THESE PEAKS EXIST AND THEY CORRELATE TO BASES, THE FOLLOWING PARAEMETERS MUST BE USED:
#### AMPLITUDE OF PEAKS, DISTANCE BETWEEN PEAKS, INTENSITY, CONFIDENCE, DERIVATIVES SIDEWAYS OF PEAK
#### MATCHES WITHIN THE ALIGNMENT CAN BE USED FOR TRANING 

# print(confidence(tmp_lst, channels_fw)[0], confidence(new_peaks, channels_fw)[0])
# print(confidence(tmp_lst, channels_fw)[1], confidence(new_peaks, channels_fw)[1])
# rename(channels_fw, guide_fw)
# print(jiggler(ploc_fw, channels_fw))
# for i in range(locs[0][0], locs[0][1]):
# for j in range(locs[1][0], locs[1][1]):
# plt.plot(list(channels_fw.values())[0], "blue")
# plt.show()