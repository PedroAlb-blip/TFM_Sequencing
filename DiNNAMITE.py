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
from Functions import filterer, peak_discovery, rename



cmd = [r"C:\Mafft\mafft.bat", "--clustalout", "--localpair", "--maxiterate", "1000", "--op", "0.8", "--ep", "2.5", "--lexp", "0.00", "--lop", "0.00", "-"]
cmd_alt = [r"C:\Mafft\mafft.bat", "--clustalout", "--localpair", "--maxiterate", "1000", "--op", "2.0", "--ep", "0.1", "--lexp", "0.00", "--lop", "0.00", "-"]
cmd_alt_alt = [r"C:\Mafft\mafft.bat", "--clustalout", "--localpair", "--maxiterate", "1000", "--op", "10.00", "--ep", "0.001", "--lexp", "0.00", "--lop", "0.00", "-"]

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
            align = f">\n{sequence_1}\n>\n{sequence_2}"
            process = subprocess.run(cmd,input=align,capture_output=True,text=True,shell=True)
            alignment_1 = str(process.stdout)
            alignment_2, alignment_3 = "", ""
            alignments = [alignment_1, alignment_2, alignment_3]
            if "*"*30 in alignment_1:
                new_file.write(" GOOD READ 1\n\n")
                new_file.write(align)
                new_file.write("\n")
                new_file.write(alignment_1)
                new_file.write("\n\n\n")
            else:
                process = subprocess.run(cmd_alt,input=align,capture_output=True,text=True,shell=True)
                alignment_2 = str(process.stdout)
                if "*"*30 in alignment_2:
                    new_file.write(" GOOD READ 2\n\n")
                    new_file.write(align)
                    new_file.write("\n")
                    new_file.write(alignment_2)
                    new_file.write("\n\n\n")
                    alignments = [alignment_1, alignment_2, alignment_3]
                else:
                    process = subprocess.run(cmd_alt_alt,input=align,capture_output=True,text=True,shell=True)
                    alignment_3 = str(process.stdout)
                    if "*"*30 in alignment_3:
                        new_file.write(" GOOD READ 3\n\n")
                        new_file.write(align)
                        new_file.write("\n")
                        new_file.write(alignment_3)
                        new_file.write("\n\n\n")
                        alignments = [alignment_1, alignment_2, alignment_3]
                    else:
                        alignments = [alignment_1, alignment_2, alignment_3]
                        aster = [aln.count("**********") for aln in alignments]
                        alignment_final = alignments[aster.index(max(aster))]
                        alignments = [alignment_1, alignment_2, alignment_3, alignment_final]
                        new_file.write(f" BAD READ {str(aster.index(max(aster)) + 1)} \n\n")
                        new_file.write(align)
                        new_file.write("\n")
                        new_file.write(alignment_final)
                        new_file.write("\n\n\n")
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
    data_file.close()

print(time.time() - start_time)
#### TO TRAIN A NN TO CHECK IF THESE PEAKS EXIST AND THEY CORRELATE TO BASES, THE FOLLOWING PARAEMETERS MUST BE USED:
#### AMPLITUDE OF PEAKS, DISTANCE BETWEEN PEAKS, INTENSITY, CONFIDENCE, DERIVATIVES SIDEWAYS OF PEAK
#### MATCHES WITHIN THE ALIGNMENT CAN BE USED FOR TRANING 