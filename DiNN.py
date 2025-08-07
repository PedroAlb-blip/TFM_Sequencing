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
import numpy as np
from Functions import filterer, peak_discovery, rename

cmd = [r"C:\Mafft\mafft.bat", "--clustalout", "--localpair", "--maxiterate", "1000", "--op", "5.0", "--ep", "2.0", "--lexp", "0.00", "--lop", "0.00", "-"]

all_files = os.listdir()

file = open("weights_2.txt", "r")
read = file.read()
arr = "".join(read.split('\n       '))

model = keras.models.load_model("c:\\Users\\Pedro\\DCYFR\\weights_checkpoint.keras")

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
        #     cumm_str_1 = ''
        #     cumm_str_2 = ''
        #     start = False
        #     end = False
        #     for j, k in zip(range(0, len(pred), 10), range(-1, -len(pred), -10)):
        #         temp_str_1 = "".join(pred[j:j+10])
        #         temp_str_2 = "".join(pred[k-10:k])
        #         if temp_str_1.count('1') > 5 and start == False:
        #             cumm_str_1 = cumm_str_1 + temp_str_1.replace('0', '1')
        #             cutoff_1 = j + 10
        #         else:
        #             start = True
        #         if temp_str_2.count('1') > 5 and end == False:
        #             cumm_str_2 = cumm_str_2 + temp_str_2.replace('0', '1')
        #             cutoff_2 = k - 10
        #         else:
        #             end = True
        #     cumm_str = cumm_str_1 + "".join(pred[cutoff_1:cutoff_2]) + cumm_str_2
        #     predic = [int(k) for k in cumm_str]
            # new_file.write(str(predic))
            # new_file.write("\n\n")
            sequence = ""
            pred = [int(num) for num in pred]
            for nt in range(0, len(pred)):
                if pred[nt] == 0:
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

#### TO TRAIN A NN TO CHECK IF THESE PEAKS EXIST AND THEY CORRELATE TO BASES, THE FOLLOWING PARAEMETERS MUST BE USED:
#### AMPLITUDE OF PEAKS, DISTANCE BETWEEN PEAKS, INTENSITY, CONFIDENCE, DERIVATIVES SIDEWAYS OF PEAK
#### MATCHES WITHIN THE ALIGNMENT CAN BE USED FOR TRANING 