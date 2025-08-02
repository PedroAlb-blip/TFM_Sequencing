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
import os
import keras
import subprocess
from keras.models import Sequential
from keras.layers import Dense, Activation
import matplotlib.pyplot as plt
import numpy as np
from Functions import filterer, peak_discovery, rename

cmd = [r"C:\Mafft\mafft.bat", "--clustalout", "--localpair", "--maxiterate", "1000", "--op", "2.0", "--ep", "0.1", "-"]

all_files = os.listdir()
file_res = open("results.txt", "r")
file_par = open("training_params.txt", "r")
lines_res = file_res.read()
lines_par = file_par.read()
results = list(filter(('').__ne__, lines_res.split('\n')))
params = list(filter(('').__ne__, lines_par.split('\n')))
ones = []
train = []
k = 0
for i in range(0, len(results)):
    if results.index(results[i]) % 3 == 0:
        ones = ones + list(ast.literal_eval(results[i]))
        train = train + list(ast.literal_eval(params[int(i/3)]))
        if len(ones) < len(train):
            ones = ones + [1]*(len(train)-len(ones))
print(len(ones), len(train))     
trainX = np.array(train)
trainY = np.array(ones)


file = open("weights_2.txt", "r")
read = file.read()
arr = "".join(read.split('\n       '))
# weight_1 = np.array(ast.literal_eval(arr.split('array(')[1].split('dtype=')[0].replace("],\n", "]").replace("\n", "")), dtype=np.float32)
# bias_1 = np.array(ast.literal_eval(arr.split('array(')[2].split('dtype=')[0].replace("],\n", "]").replace("\n", "")), dtype=np.float32)
# weight_2 = np.array(ast.literal_eval(arr.split('array(')[3].split('dtype=')[0].replace("],\n", "]").replace("\n", "")), dtype=np.float32)
# bias_2 = np.array(ast.literal_eval(arr.split('array(')[4].split('dtype=')[0].replace("],\n", "]").replace("\n", "")), dtype=np.float32)
# weight_3 = np.array(ast.literal_eval(arr.split('array(')[5].split('dtype=')[0].replace("],\n", "]").replace("\n", "")), dtype=np.float32)
# bias_3 = np.array(ast.literal_eval(arr.split('array(')[6].split('dtype=')[0].replace("],\n", "]").replace("\n", "")), dtype=np.float32)

# model = Sequential()
# model.add(Dense(32, input_dim=15)) 
# model.add(Activation('relu'))
# model.add(Dense(16)) 
# model.add(Activation('relu'))
# model.add(Dense(1))
# model.add(Activation('sigmoid'))

model = keras.models.load_model("c:\\Users\\Pedro\\DCYFR\\weights_checkpoint.keras")

# j = 1
# for layer in [model.layers[0], model.layers[2], model.layers[4]]:
#     layer.set_weights([eval(f"weight_{j}")[0], eval(f"bias_{j}")[0]])
#     j = j + 1

# checkpoint_path = "c:\\Users\\Pedro\\DCYFR\\weights_checkpoint.keras"
# checkpoint = keras.callbacks.ModelCheckpoint(filepath=checkpoint_path, save_freq="epoch", monitor='accuracy', mode='max', save_best_only=True)
# model.compile(loss=keras.losses.BinaryFocalCrossentropy(gamma=2.0, alpha=0.25), optimizer=keras.optimizers.Adamax(learning_rate=0.001), metrics=['accuracy']) #
# history = model.fit(trainX, trainY, epochs=100, batch_size=80, verbose=1, validation_split=0.25, callbacks=[checkpoint])

# keras.models.load_model(checkpoint_path)

# weighty_file = open("weights_3.txt", "w")
# weighty_file.write(str(model.get_weights()))
# plt.plot(history.history['accuracy'])
# plt.plot(history.history['val_accuracy'])
# plt.title('model accuracy')
# plt.ylabel('accuracy')
# plt.xlabel('epoch')
# plt.legend(['train', 'test'], loc='upper left')
# plt.show()
# plt.plot(history.history['loss'])
# plt.plot(history.history['val_loss'])
# plt.title('model loss')
# plt.ylabel('loss')
# plt.xlabel('epoch')
# plt.legend(['train', 'test'], loc='upper left')
# plt.show()

new_file = open("attempt.txt", "a")
for i in all_files:
    if i.startswith("Sp"):
        data_file = open(i, "r")
        new_file.write(str(i.split(".txt")[0]))
        new_file.write("\n\n")
        file_read = data_file.read()
        every = list(filter(('').__ne__, file_read.split('\n')))
        channels_fw, channels_rv, guide_fw, guide_rv = ast.literal_eval(every[10]), ast.literal_eval(every[12]), every[14], every[16]
        channels_fw = rename(channels_fw, guide_fw)
        channels_rv = rename(channels_rv, guide_rv)
        is_fw = np.array(filterer(channels_fw, peak_discovery(channels_fw)[0]))
        is_rv = np.array(filterer(channels_rv, peak_discovery(channels_rv)[0]))
        res_fw = np.ndarray.tolist(model.predict(is_fw))
        res_rv = np.ndarray.tolist(model.predict(is_rv))
        seq_fw = peak_discovery(channels_fw)[2]
        seq_rv = peak_discovery(channels_rv)[2]
        pred_fw = [str(int(round(v[0], 0))) for v in res_fw[:]]
        pred_rv = [str(int(round(v[0], 0))) for v in res_rv[:]]
        sequence_1, sequence_2 = "", ""
        for pred, seq in zip((pred_fw, pred_rv), (seq_fw, seq_rv)):
            cumm_str_1 = ''
            cumm_str_2 = ''
            start = False
            end = False
            for j, k in zip(range(0, len(pred), 10), range(-1, -len(pred), -10)):
                temp_str_1 = "".join(pred[j:j+10])
                temp_str_2 = "".join(pred[k-10:k])
                if temp_str_1.count('1') > 5 and start == False:
                    cumm_str_1 = cumm_str_1 + temp_str_1.replace('0', '1')
                    cutoff_1 = j + 10
                else:
                    start = True
                if temp_str_2.count('1') > 5 and end == False:
                    cumm_str_2 = cumm_str_2 + temp_str_2.replace('0', '1')
                    cutoff_2 = k - 10
                else:
                    end = True
            cumm_str = cumm_str_1 + "".join(pred[cutoff_1:cutoff_2]) + cumm_str_2
            predic = [int(k) for k in cumm_str]
            # new_file.write(str(predic))
            # new_file.write("\n\n")
            sequence = ""
            for nt in range(0, len(predic)):
                if predic[nt] == 0:
                    sequence = sequence + seq[nt]
            # new_file.write(str(sequence))
            # new_file.write("\n\n")
            if sequence_1 == "":
                sequence_1 = sequence
            else:
                sequence_2 = sequence
        align = f">\n{sequence_1}\n>\n{sequence_2}"
        process = subprocess.run(cmd,input=align,capture_output=True,text=True,shell=True)
        align = process.stdout
        new_file.write(str(align))
        new_file.write("\n\n\n")
        data_file.close()
# file_1 = open("Sp113-VP7a.txt", "r")
# file_2 = open("Sp116-VP7a.txt", "r")
# data = file_1.read()
# datt = file_2.read()
# all = list(filter(('').__ne__, data.split('\n')))
# ali = list(filter(('').__ne__, datt.split('\n')))
# c_f_13 = ast.literal_eval(all[10])
# c_r_13 = ast.literal_eval(all[12])
# g_f = all[14]
# g_r = all[16]
# c_f_16 = ast.literal_eval(ali[10])
# c_r_16 = ast.literal_eval(ali[12])

# c_f_13 = rename(c_f_13, g_f)
# c_r_13 = rename(c_r_13, g_r)
# c_f_16 = rename(c_f_16, g_f)
# c_r_16 = rename(c_r_16, g_r)

# dict = {"c_f_13" : c_f_13, "c_r_13" : c_r_13, "c_f_16" : c_f_16, "c_r_16" : c_r_16}
# keys = list(dict.keys())
# for i in range(0, 4):
#     dataPred = model.predict(np.array(filterer(dict[keys[i]], peak_discovery(dict[keys[i]])[0])))
#     pred = [int(round(k, 0)) for k in list(np.ndarray.flatten(dataPred))]
#     print(pred, '\n\n')

#### TO TRAIN A NN TO CHECK IF THESE PEAKS EXIST AND THEY CORRELATE TO BASES, THE FOLLOWING PARAEMETERS MUST BE USED:
#### AMPLITUDE OF PEAKS, DISTANCE BETWEEN PEAKS, INTENSITY, CONFIDENCE, DERIVATIVES SIDEWAYS OF PEAK
#### MATCHES WITHIN THE ALIGNMENT CAN BE USED FOR TRANING 