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

file_res = open("results.txt", "r")
file_par = open("training_params.txt", "r")
lines_res = file_res.read()
lines_par = file_par.read()
results = list(filter(('').__ne__, lines_res.split('\n')))
params = list(filter(('').__ne__, lines_par.split('\n')))
ones = []
train = []
for i in range(0, len(results)):
    if results.index(results[i]) % 3 == 2:
            ones = ones + list(ast.literal_eval(results[i]))
            train = train + list(ast.literal_eval(params[int(round(i/3,0)) - 1]))
            if len(ones) < len(train):
                ones = ones + [1]*(len(train)-len(ones))
trainX = np.array(train)
trainY = np.array(ones)

model = Sequential()
file = open("weights.txt", "r")
read = file.read()
arr = "".join(read.split('\n       '))
array_1 = np.array(ast.literal_eval(arr.split('array(')[1].split('dtype=')[0].replace("],\n", "]").replace("\n", "")), dtype=np.float32)
array_2 = np.array(ast.literal_eval(arr.split('array(')[2].split('dtype=')[0].replace("],\n", "]").replace("\n", "")), dtype=np.float32)
array_3 = np.array(ast.literal_eval(arr.split('array(')[3].split('dtype=')[0].replace("],\n", "]").replace("\n", "")), dtype=np.float32)
array_4 = np.array(ast.literal_eval(arr.split('array(')[4].split('dtype=')[0].replace("],\n", "]").replace("\n", "")), dtype=np.float32)
array_5 = np.array(ast.literal_eval(arr.split('array(')[5].split('dtype=')[0].replace("],\n", "]").replace("\n", "")), dtype=np.float32)
array_6 = np.array(ast.literal_eval(arr.split('array(')[6].split('dtype=')[0].replace("],\n", "]").replace("\n", "")), dtype=np.float32)
print(array_5)
model.add(Dense(20, input_dim=10, activation='relu', weights = [array_1, array_2]))
model.add(Dense(15, activation='sigmoid', weights = [array_3, array_4]))
model.add(Dense(1, activation='sigmoid', weights = [array_5[0], array_6[0]]))
# model.compile(loss='mean_squared_error', optimizer='adam')
# model.fit(trainX, trainY, epochs=200, batch_size=32, verbose=1)
# weights.write(str(model.get_weights()))

file_1 = open("Sp113-VP7a.txt", "r")
file_2 = open("Sp116-VP7a.txt", "r")
data = file_1.read()
datt = file_2.read()
all = list(filter(('').__ne__, data.split('\n')))
ali = list(filter(('').__ne__, datt.split('\n')))
c_f_13 = ast.literal_eval(all[10])
c_r_13 = ast.literal_eval(all[12])
g_f = all[14]
g_r = all[16]
c_f_16 = ast.literal_eval(ali[10])
c_r_16 = ast.literal_eval(ali[12])

c_f_13 = rename(c_f_13, g_f)
c_r_13 = rename(c_r_13, g_r)
c_f_16 = rename(c_f_16, g_f)
c_r_16 = rename(c_r_16, g_r)

dict = {"c_f_13" : c_f_13, "c_r_13" : c_r_13, "c_f_16" : c_f_16, "c_r_16" : c_r_16}
keys = list(dict.keys())
for i in range(0, 4):
    dataPred = model.predict(np.array(filterer(dict[keys[i]], peak_discovery(dict[keys[i]])[0])[0]))
    pred = list(np.ndarray.flatten(dataPred))
    print(pred, '\n\n')
# dataPrediction = model.predict(np.array(filterer(channels_fw, peak_discovery(channels_fw)[0])))
# pred = list(np.ndarray.flatten(dataPrediction))
# print(pred)
# print('\n\n')



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