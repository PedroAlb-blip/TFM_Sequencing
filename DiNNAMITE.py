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

new_file = open("alignments.txt", "a")
other_new_file = open("consensi.txt", "a")

for i in all_files: #Iterating over all files in the sequences folder
    data_file = open(f"Sequences\\{i}", "r")
    new_file.write(str(i.split(".txt")[0]))
    other_new_file.write(str(i.split(".txt")[0])) # Adds the sequence name to the output files
    file_read = data_file.read()
    every = list(filter(('').__ne__, file_read.split('\n')))
    channels_fw, channels_rv, guide_fw, guide_rv = ast.literal_eval(every[10]), ast.literal_eval(every[12]), every[14], every[16]
    channels_fw = rename(channels_fw, guide_fw)
    channels_rv = rename(channels_rv, guide_rv)
    correction = 0
    while correction < 3:
        sequence_1, sequence_2 = "", ""
        lod, ploc_sub = [], []
        for channels in (channels_fw, channels_rv): ### The sequences are clipped eliminitaing untruthful peaks i.e. 1s in the binary vector
            is_seq = np.array(filterer(channels, peak_discovery(channels)[0]))
            res_seq = np.ndarray.tolist(model.predict(is_seq)) ### Here the neural network is called to produce a binary vector
            pred = [str(int(round(v[0], 0))) for v in res_seq[:]] ### Here the values are rounded to either 1 or 0 by mathematical rounding
            # new_file.write(str(pred))
            # new_file.write("\n\n")
            seq = peak_discovery(channels)[2]
            sequence = ""
            pred = [int(num) for num in pred]
            ploc = peak_discovery(channels)[1] ### The full list is kept
            true_ploc = []
            for nt in range(0, len(pred)):
                if pred[nt] == 0:
                    sequence = sequence + seq[nt] ### the sequence is generated
                    true_ploc.append(ploc[nt]) ### True peaks are added to a list
            # new_file.write(str(sequence))
            # new_file.write("\n\n")
            if sequence_1 == "":
                sequence_1 = sequence
                channels_1 = channels
                ploc_1 = true_ploc ### Sequence and channels are added to either 1 for forward or 2 for reverse
            elif sequence_1 != "" and sequence_2 == "":
                sequence_2 = sequence
                channels_2 = channels
                ploc_2 = true_ploc
            if sequence_1 != "" and sequence_2 != "":
                alignment = ""
                alignments = []
                cover = 1
                delay_fw, delay_rv = 1, 1
                if len(sequence_1) > len(sequence_2):
                    delay_fw = len(sequence_2)/len(sequence_1) 
                else:
                    delay_rv = len(sequence_1)/len(sequence_2)
                while "*"*25 not in alignment and cover > 0.1: ### 25 consecutive * serve as a marker for good alignment and the minimal coverage of secuence is set to 10%
                    align = f">\n{sequence_1[int(round(((len(sequence_1)-1)-(len(sequence_1)-1)*cover)*delay_fw, 0)) : -1]}\n>\n{sequence_2[0:int(round(((len(sequence_2)-1)*cover)*delay_rv, 0))]}"
                    process = subprocess.run(cmd,input=align,capture_output=True,text=True,shell=True)
                    alignment = str(process.stdout)
                    alignments.append(alignment)
                    cover = cover/1.25 ### As the tool finds no good alignment it reduces the amount of sequence covered
                choice = [aln.count("*****") for aln in alignments] ### Of all sequences the one with the most accounts of 5 consecutive matches are taken to be the better ones
                if max(choice) == 0: ### In case no good alignment is found then the last is chosen, with the least coverage
                    align_split = alignments[-1].split('\n')
                    max_choice = len(choice) - 1
                else:
                    align_split = alignments[choice.index(max(choice))].split('\n')
                    max_choice = choice.index(max(choice))
                align_purged = [u for u in align_split if len(u) > 1]
                new_file.write(f" NUMBER OF READ: {len(alignments)}; COVERAGE {int(cover*(1.25**(len(choice) - max_choice))*100)}% ; CORRECTION LEVEL {correction}\n\n")
                asterisk, seq_fw_al, seq_rv_al = "", "", ""
                for chunk in range(0,len(align_purged)): ### Alignment is chopped into 3 parts: forward, reverse and alignment and they are filled
                    if chunk % 3 == 0 and chunk != 0:
                        asterisk = asterisk + align_purged[chunk][16:].replace("."," ")
                    elif chunk % 3 == 1:
                        seq_fw_al = seq_fw_al + align_purged[chunk][16:]
                    elif chunk % 3 == 2:
                        seq_rv_al = seq_rv_al + align_purged[chunk][16:]
                ### Sequences are completed with the parts not present in the alignment
                seq_fw_al_ad = sequence_1[0:int(round((((len(sequence_1)-1)-(len(sequence_1)-1)*cover*1.25**(len(choice) - choice.index(max(choice))))*delay_fw),0))]
                seq_rv_al_ad = sequence_2[int(round(((len(sequence_2)-1)*cover*1.25**(len(choice) - choice.index(max(choice))))*delay_rv, 0)): -1]
                asterisk_1 = "-"*len(seq_fw_al_ad)
                asterisk_2 = "-"*len(seq_rv_al_ad)
                larger = [len(splitted) for splitted in asterisk.split(" ")] ### List containing all asterisk regions
                asterisk = asterisk_1 + asterisk + asterisk_2
                seq_fw_al = seq_fw_al_ad + seq_fw_al + asterisk_2
                seq_rv_al = asterisk_1 + seq_rv_al + seq_rv_al_ad
                align = f">FORWARD\n{sequence_1}\n>REVERSE\n{sequence_2}\n\n"
                new_file.write(align)
                new_file.write("\n")
                ### The full alignments are printed into the files
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
                ### Here the correction begins, then it can be iterated over the sequence 1 and sequence 2 if statement, so that the alignments are printed
                gaps = re.finditer(r"\*{5} {1,10}\*{3,}", asterisk) ### The search for gaps is only performed in non stochastic matches
                indices = [[match.start(), match.end()] for match in gaps]
                for index in indices:
                    channels_proxy_fw, channels_proxy_rv = {key:[] for key in channels_fw.keys()}, {key:[] for key in channels_fw.keys()}
                    fw_index_1 = ploc_1[sequence_1.find(seq_fw_al[index[0]:index[1]].replace("-", "").upper())]
                    fw_index_2 = ploc_1[sequence_1.find(seq_fw_al[index[0]:index[1]].replace("-", "").upper()) + index[1] - index[0] - seq_fw_al[index[0]:index[1]].count("-")]
                    rv_index_1 = ploc_2[sequence_2.find(seq_rv_al[index[0]:index[1]].replace("-", "").upper())]
                    rv_index_2 = ploc_2[sequence_2.find(seq_rv_al[index[0]:index[1]].replace("-", "").upper()) + index[1] - index[0] - seq_rv_al[index[0]:index[1]].count("-")]
                    offset_fw = (2**((len(asterisk) - index[0])/len(asterisk)))/((2**(index[0]/len(asterisk)))+2**((len(asterisk) - index[0])/len(asterisk)))
                    offset_rv = 2**(index[0]/len(asterisk))/((2**(index[0]/len(asterisk)))+2**((len(asterisk) - index[0])/len(asterisk)))
                    for key in channels_fw.keys():
                        channels_proxy_fw[key] = list(channels_fw[key][fw_index_1:fw_index_2])
                        channels_proxy_rv[key] = list(channels_rv[key][rv_index_1:rv_index_2])
                    fw_values = list(channels_proxy_fw.values())
                    rv_values = list(channels_proxy_rv.values())
                    channels_fw_sum = [x + y + z + t for x, y, z, t in zip(fw_values[0], fw_values[1], fw_values[2], fw_values[3])]
                    channels_rv_sum = [x + y + z + t for x, y, z, t in zip(rv_values[0], rv_values[1], rv_values[2], rv_values[3])]
                    if fw_index_2 - fw_index_1 > rv_index_2 - rv_index_1:
                        length = index[1] - index[0] - seq_fw_al[index[0]:index[1]].count("-")
                        step = 0
                        while step < fw_index_2 - fw_index_1 - (rv_index_2 - rv_index_1):
                            roll = step % length
                            channels_slice = channels_fw_sum[int(round(roll*(len(channels_fw_sum))/length, 0)):int(round((roll+1)*(len(channels_fw_sum))/length, 0))]
                            bye = channels_slice.index(min(channels_slice)) + int(round(roll*(len(channels_fw_sum))/length,0))
                            channels_fw_sum.pop(bye)
                            for key in channels_proxy_fw.keys():
                                try:
                                    channels_proxy_fw[key][bye+1] = (channels_proxy_fw[key][bye] + channels_proxy_fw[key][bye + 1])/2
                                except IndexError:
                                    channels_proxy_fw[key][bye-1] = (channels_proxy_fw[key][bye] + channels_proxy_fw[key][bye - 1])/2
                                channels_proxy_fw[key].pop(bye)
                            step += 1
                    elif  rv_index_2 - rv_index_1 > fw_index_2 - fw_index_1:
                        length = index[1] - index[0] - seq_rv_al[index[0]:index[1]].count("-")
                        step = 0
                        while step < rv_index_2 - rv_index_1 - (fw_index_2 - fw_index_1):
                            roll = step % length
                            channels_slice = channels_rv_sum[int(round(roll*(len(channels_rv_sum))/length, 0)):int(round((roll+1)*(len(channels_rv_sum))/length, 0))]
                            bye = channels_slice.index(min(channels_slice)) + int(round(roll*(len(channels_rv_sum))/length,0))
                            channels_rv_sum.pop(bye)
                            for key in channels_proxy_fw.keys():
                                try:
                                    channels_proxy_rv[key][bye + 1] = (channels_proxy_rv[key][bye] + channels_proxy_rv[key][bye + 1])/2
                                except IndexError:
                                    channels_proxy_rv[key][bye - 1] = (channels_proxy_rv[key][bye] + channels_proxy_rv[key][bye - 1])/2
                                channels_proxy_rv[key].pop(bye)
                            step += 1
                    else:
                        continue
                    channels_sum = {key:[] for key in channels_fw.keys()}
                    # color = ["blue", "red", "green", "orange"]
                    # p=0
                    for key in list(channels_sum.keys()):
                        channels_sum[key] = [x*offset_fw + y*offset_rv for x, y in zip(channels_proxy_fw[key], channels_proxy_rv[key])]
                    #     plt.plot(channels_sum[key], color = color[p], label = key)
                    #     channels_sum[key] = [0]*8000 + channels_sum[key] + [0]*6000
                    #     p+=1
                    # for peak in peak_discovery(channels_sum)[1]:
                    #     plt.text(peak, max([channels_sum[key][peak] for key in channels_sum.keys()]), [key for key in channels_sum.keys() if peak in list(peak_discovery(channels_sum)[0][key])][0])
                    # plt.legend(loc="upper left")
                    # plt.show()
                    lod.append(channels_sum)
                    ploc_sub.append([fw_index_1, fw_index_2, rv_index_1, rv_index_2])
        lod = lod[::-1] ### These have been flipped to avoid index shennanigans as the substitutions can be smaller
        ploc_sub = ploc_sub[::-1]
        for sub in range(0,len(lod)):
            for key in lod[sub].keys():
                channels_fw[key] = channels_fw[key][0:ploc_sub[sub][0]] + lod[sub][key] + channels_fw[key][ploc_sub[sub][1]:]
                channels_rv[key] = channels_rv[key][0:ploc_sub[sub][2]] + lod[sub][key] + channels_rv[key][ploc_sub[sub][3]:]
        correction += 1
    pos_1, pos_2 = sequence_1.find(seq_fw_al[asterisk.find("*"*max(larger)):].replace("-","").upper()), sequence_2.find(seq_rv_al[asterisk.find("*"*max(larger)):].replace("-","").upper())
    dict_1, dict_2 = {key:[] for key in channels_1.keys()}, {key:[] for key in channels_2.keys()}
    ### DOP are created within the better read regions 
    for app in range(pos_1, pos_1 + max(larger)):
        dict_1[sequence_1[app]].append(ploc_1[app])
    for app in range(pos_2, pos_2 + max(larger)):
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
    new_file.write(">> DINNAMITE consensus \n")
    new_file.write(consensus)
    other_new_file.write("\n\n>> DINNAMITE consensus \n")
    other_new_file.write(consensus)
    new_file.write("\n\n\n")
    other_new_file.write("\n\n\n")
    data_file.close()
new_file.close()
other_new_file.close()

print(time.time() - start_time)
#### TO TRAIN A NN TO CHECK IF THESE PEAKS EXIST AND THEY CORRELATE TO BASES, THE FOLLOWING PARAEMETERS MUST BE USED:
#### AMPLITUDE OF PEAKS, DISTANCE BETWEEN PEAKS, INTENSITY, CONFIDENCE, DERIVATIVES SIDEWAYS OF PEAK
#### MATCHES WITHIN THE ALIGNMENT CAN BE USED FOR TRANING 

        # alignment = [aln for aln in alignments if aln != ""][-1]
        # align_split = alignment.split('\n')
        # align_purged = [u for u in align_split if len(u) > 1]
        # only, seq_fw_al, seq_rv_al = "", "", ""
        # for chunk in range(0,len(align_purged)):
        #     if chunk % 3 == 0 and chunk != 0:
        #         only = only + align_purged[chunk][16:]
        #     elif chunk % 3 == 1:
        #         seq_fw_al = seq_fw_al + align_purged[chunk][16:]
        #     elif chunk % 3 == 2:
        #         seq_rv_al = seq_rv_al + align_purged[chunk][16:]
        # asterisk = only[re.search(r"[actg]{9}", seq_fw_al).start() : - re.search(r"[actg]{9}", seq_rv_al[::-1]).start() - 1]
        # seq_fw_pg = seq_fw_al[re.search(r"[actg]{9}", seq_fw_al).start() : - re.search(r"[actg]{9}", seq_rv_al[::-1]).start() - 1]
        # seq_rv_pg = seq_rv_al[re.search(r"[actg]{9}", seq_fw_al).start() : - re.search(r"[actg]{9}", seq_rv_al[::-1]).start() - 1]
        # gaps = re.finditer(r"[\*\.]{2,6} +[\*\.]{2,6}", asterisk)
        # colors = {"A" : "green", "G" : "blue", "T" : "black", "C" : "red"}
        # for gap in gaps:
        #     dict_sum = {key : [] for key in list(channels_1.keys())}
        #     for key in dict_sum:
        #         pos_1 = 2**((len(asterisk)-((gap.start() + gap.end())/2))/len(asterisk))
        #         pos_2 = 2**(((gap.start() + gap.end())/2)/len(asterisk))
        #         ratio_1 = pos_1/(pos_1 + pos_2)
        #         ratio_2 = pos_2/(pos_1 + pos_2)
        #         list_1 = channels_1[key][ploc_1[sequence_1.find(seq_fw_pg[gap.start() - 4 :gap.end() + 4].replace("-", "").upper())] : ploc_1[sequence_1.find(seq_fw_pg[gap.start() - 4 :gap.end() + 4].replace("-", "").upper()) + len(seq_fw_pg[gap.start() - 4 :gap.end() + 4].replace("-", "")) - 1]]
        #         list_2 = channels_2[key][ploc_2[sequence_2.find(seq_rv_pg[gap.start() - len(seq_rv_pg) - 4 : gap.end() - len(seq_rv_pg) + 4].replace("-","").upper())] : ploc_2[sequence_2.find(seq_rv_pg[gap.start() - len(seq_rv_pg) - 4 : gap.end() - len(seq_rv_pg) + 4].replace("-","").upper()) + len(seq_rv_pg[gap.start() - len(seq_rv_pg) - 4 : gap.end() - len(seq_rv_pg) + 4].replace("-","")) - 1]]
        #         dict_sum[key] = tuple([0]*800) + tuple([(list_1[ind]*ratio_1 + list_2[ind]*ratio_2)/2 for ind in range(0, min(len(list_1), len(list_2)))]) + tuple([0]*2000)
        #         plt.plot(dict_sum[key], color = colors[key])
        #     is_gap = np.array(filterer(dict_sum, peak_discovery(dict_sum)[0]))
        #     seq_gap = peak_discovery(dict_sum)[2]
        #     gap_pred = np.ndarray.tolist(model.predict(is_gap))
        #     gap_bin = [int(round(num[0], 0)) for num in gap_pred[:]]
        #     print(gap_bin)
        #     temp_seq = ''
        #     for bin in range(0, len(gap_bin)):
        #         if gap_bin[bin] == int(0):
        #             temp_seq = temp_seq + seq_gap[bin]
        #     print(seq_fw_pg[gap.start() - 4 :gap.end() + 4].replace("-", "").upper(), "\n", seq_rv_pg[gap.start() - len(seq_rv_pg) - 4 : gap.end() - len(seq_rv_pg) + 4].replace("-", "").upper(), "\n", temp_seq)
        #     plt.show()
