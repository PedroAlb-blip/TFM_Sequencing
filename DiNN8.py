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
import time
start_time = time.time()
import ast
import matplotlib.pyplot as plt
import os
import keras
import subprocess
import numpy as np
import re
from Functions import filterer, peak_discovery, rename, confidence



cmd = [r"C:\Mafft\mafft.bat", "--clustalout", "--localpair", "--maxiterate", "500", "--op", "0.8", "--ep", "2.5", "--lexp", "6.00", "--lop", "6.00", "-"]
### ADD A CONDITION TO REPEAT UNDER DIFFERENT PARAMETERS IF THE ALIGNMENT LOOKS LIKE SHIT FOR INSTANCE LESS THAN 20 ASTERISKS TOGETHER

all_files = os.listdir("c:\\Users\\Pedro\\DCYFR\\Sequences")

model = keras.models.load_model("c:\\Users\\Pedro\\DCYFR\\weights_checkpoint.keras")

new_file = open("attempt.txt", "a")
for i in all_files:
    data_file = open(f"Sequences\\{i}", "r")
    new_file.write(str(i.split(".txt")[0]))
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
    for pred, seq, channels in zip((pred_fw, pred_rv), (seq_fw, seq_rv), (channels_fw, channels_rv)):
        # new_file.write(str(predic))
        # new_file.write("\n\n")
        sequence = ""
        pred = [int(num) for num in pred]
        ploc = peak_discovery(channels)[1]
        true_ploc = []
        for nt in range(0, len(pred)):
            if pred[nt] == 0:
                sequence = sequence + seq[nt]
                true_ploc.append(ploc[nt])
        # new_file.write(str(sequence))
        # new_file.write("\n\n")
        if sequence_1 == "":
            sequence_1 = sequence
            channels_1 = channels
            ploc_1 = true_ploc
        else:
            sequence_2 = sequence
            channels_2 = channels
            ploc_2 = true_ploc
            alignment = ""
            alignments = []
            cover = 1
            offset_fw, offset_rv = 1, 1
            if len(sequence_1) > len(sequence_2):
                offset_fw = len(sequence_2)/len(sequence_1) 
            else:
                offset_rv = len(sequence_1)/len(sequence_2)
            while "*"*25 not in alignment and cover > 0.1:
                align = f">\n{sequence_1[int(round(((len(sequence_1)-1)-(len(sequence_1)-1)*cover)*offset_fw, 0)) : -1]}\n>\n{sequence_2[0:int(round(((len(sequence_2)-1)*cover)*offset_rv, 0))]}"
                process = subprocess.run(cmd,input=align,capture_output=True,text=True,shell=True)
                alignment = str(process.stdout)
                alignments.append(alignment)
                cover = cover/1.25
            choice = [aln.count("*****") for aln in alignments]
            align_split = alignments[choice.index(max(choice))].split('\n')
            align_purged = [u for u in align_split if len(u) > 1]
            new_file.write(f" NUMBER OF READ: {len(alignments)}; COVERAGE {int(cover*(len(choice) - choice.index(max(choice)))*125)}% \n\n")
            asterisk, seq_fw_al, seq_rv_al = "", "", ""
            for chunk in range(0,len(align_purged)):
                if chunk % 3 == 0 and chunk != 0:
                    asterisk = asterisk + align_purged[chunk][16:].replace("."," ")
                elif chunk % 3 == 1:
                    seq_fw_al = seq_fw_al + align_purged[chunk][16:]
                elif chunk % 3 == 2:
                    seq_rv_al = seq_rv_al + align_purged[chunk][16:]
            seq_fw_al_ad = sequence_1[0:int(round((((len(sequence_1)-1)-(len(sequence_1)-1)*cover*(len(choice) - choice.index(max(choice)))*1.25)*offset_fw),0))]
            seq_rv_al_ad = sequence_2[int(round(((len(sequence_2)-1)*(len(choice) - choice.index(max(choice)))*cover*1.25)*offset_rv, 0)): -1]
            asterisk_1 = "-"*len(seq_fw_al_ad)
            asterisk_2 = "-"*len(seq_rv_al_ad)
            larger = [len(splitted) for splitted in asterisk.split(" ")]
            asterisk = asterisk_1 + asterisk + asterisk_2
            seq_fw_al = seq_fw_al_ad + seq_fw_al + asterisk_2
            seq_rv_al = asterisk_1 + seq_rv_al + seq_rv_al_ad
            align = f">FORWARD\n{sequence_1}\n>REVERSE\n{sequence_2}\n\n"
            new_file.write(align)
            new_file.write("\n")
            for line in range(0,len(asterisk),60):
                if line + 60 > len(asterisk):
                    new_file.write(seq_fw_al[line:])
                    new_file.write("\n")
                    new_file.write(seq_rv_al[line:])
                    new_file.write("\n")
                    new_file.write(asterisk[line:])
                    new_file.write("\n\n")
                    break
                new_file.write(seq_fw_al[line:line+60])
                new_file.write("\n")
                new_file.write(seq_rv_al[line:line+60])
                new_file.write("\n")
                new_file.write(asterisk[line:line+60])
                new_file.write("\n\n")
            pos_1, pos_2 = sequence_1.find(seq_fw_al[asterisk.find("*"*max(larger)):].replace("-","").upper()), sequence_2.find(seq_rv_al[asterisk.find("*"*max(larger)):].replace("-","").upper())
            dict_1, dict_2 = {key:[] for key in channels_1.keys()}, {key:[] for key in channels_2.keys()}
            for app in range(pos_1, pos_1+len(seq_fw_al[asterisk.find("*"*max(larger)):len(sequence_1)-1])):
                dict_1[sequence_1[app]].append(ploc_1[app])
            for app in range(pos_2, pos_2+len(seq_rv_al[asterisk.find("*"*max(larger)):len(sequence_2)-1])):
                dict_2[sequence_2[app]].append(ploc_2[app])
            conf_1, conf_2 = confidence(dict_1, channels_1), confidence(dict_2, channels_2)
            breaker = 0
            for confis in range(0,min(len(conf_1[0]), len(conf_2[0]))):
                ratio = (confis)/min(len(conf_1[0]), len(conf_2[0]))
                switch_1 = conf_1[0][confis] * conf_1[1][confis]
                switch_2 = conf_2[0][confis] * conf_2[1][confis]
                if switch_2 > switch_1 and breaker != 5:
                    breaker = breaker + 1
                elif switch_1 > switch_2:
                    breaker = 0
                elif breaker == 5 or confis == min(len(conf_1[0]), len(conf_2[0])) - 1:
                    pos_1 = pos_1 + confis - 5
                    pos_2 = pos_2 + confis - 5
                    break
            consensus = sequence_1[:pos_1] + sequence_2[pos_2:]
            new_file.write(">> DINN8ing consensus \n")
            new_file.write(consensus)
            new_file.write("\n\n\n")
    data_file.close()

print(time.time() - start_time)
#### TO TRAIN A NN TO CHECK IF THESE PEAKS EXIST AND THEY CORRELATE TO BASES, THE FOLLOWING PARAEMETERS MUST BE USED:
#### AMPLITUDE OF PEAKS, DISTANCE BETWEEN PEAKS, INTENSITY, CONFIDENCE, DERIVATIVES SIDEWAYS OF PEAK
#### MATCHES WITHIN THE ALIGNMENT CAN BE USED FOR TRANING 