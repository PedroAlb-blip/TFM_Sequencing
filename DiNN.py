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
from Functions import filterer, peak_discovery, rename



#### True is 0 and False is 1
file = open("Sp101b-VP7.txt", "r")
data = file.read()
all = list(filter(('').__ne__, data.split('\n')))
###### IMPORTANT INFORMATION IN: 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22 and 24
channels_fw = ast.literal_eval(all[10])
channels_rv = ast.literal_eval(all[12])
guide_fw = all[14]
guide_rv = all [16]
channels_fw = rename(channels_fw, guide_fw)
channels_rv = rename(channels_rv, guide_rv)
locs = ast.literal_eval(all[24])
align = all[18]
fw_seq = all[2]
rv_seq = all[4]


train = filterer(channels_fw, peak_discovery(channels_fw))
print(train)
trainX = np.array(train[0])
trainY = np.array(train[1])
model = Sequential()

model.add(Dense(8, input_dim=8, activation='relu'))
model.add(Dense(1, activation='sigmoid'))
model.compile(loss='mean_squared_error', optimizer='adam')
model.fit(trainX, trainY, epochs=200, batch_size=32, verbose=1)

dataPrediction = model.predict(np.array(filterer(channels_fw, peak_discovery(channels_fw))))
pred = list(np.ndarray.flatten(dataPrediction))
print(pred)
print('\n\n')
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