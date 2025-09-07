import os
import subprocess
from Bio import SeqIO
import copy

#### PERHAPS A CONFIDENCE METRIC CAN BE CALCULATED AS THE VALUE OF EACH CHANNEL DIVIDED BY THE TOTAL AND CHOOSE ONWARDS FROM A THRESHOLD OF HIGH REPETITIONS
def confidence(dop, dol): ### DOP is a dictionary of peaks that is the peaks in their respective channels, DOL is a dictionary of lists that is the whole channel lists
    conl = {key : [] for key in list(dop.keys())} 
    intens = {key : [] for key in list(dop.keys())}
    conf, ints, peaks, peak_dist_fw, peak_dist_bw = [], [], [], [], []
    for keys in list(dop.keys()):
        for i in dop[keys]:
            peaks.append(i) ### appended with the peak position
            sums = sum([dol[c][i] for c in dol.keys()]) ### appended with the sum of intensities at peak position
            intens[keys].append(dol[keys][i]) ### appended with all intensities at peak position
            try:
                conme = dol[keys][i]/sums ### confidence is calculated as the intensity of each channel over the whole
            except ZeroDivisionError:
                conme = 0
            conl[keys].append(conme)
    peaks.sort()
    for i in range(0, len(peaks)):
        ### Here peak distance is calculated to the next and previous peaks
        try:
            peak_dist_fw.append(peaks[i + 1] - peaks [i])
        except IndexError:
            peak_dist_fw.append(-1)
        try:
            peak_dist_bw.append(peaks [i] - peaks[i - 1]) 
        except i == 0:
            peak_dist_bw.append(-1)
        ### Here another neat trick is used, to avoid possible .index() mistakes the values are eliminated from the lists as they are transcribed to the finalised lists
        try:
            vals = [dop[key][0] for key in list(dop.keys())]
        except IndexError: ### As we are eliminating values within lists a point will come of error at dop[key][0] as it does not exist
            leng = [len(dop[key]) for key in list(dop.keys())]
            key = list(dop.keys())[leng.index(min(leng))]
            dop.pop(key, None)
            conl.pop(key, None)
            intens.pop(key, None)
            vals = [dop[key][0] for key in list(dop.keys())]
        conf.append(conl[list(dop.keys())[vals.index(min(vals))]][0])
        ints.append(intens[list(dop.keys())[vals.index(min(vals))]][0])
        conl[list(dop.keys())[vals.index(min(vals))]].pop(0)
        intens[list(dop.keys())[vals.index(min(vals))]].pop(0)
        dop[list(dop.keys())[vals.index(min(vals))]].pop(0)
    return conf, ints, peak_dist_fw, peak_dist_bw

def rename(dol, guide):
    ### Easy function here, it renames the channels (DOL) to the corresponding nucleotides as stated in the guides
    tmp_lst = list(dol.keys())
    for i in range(0, len(dol.keys())):
        dol[guide[i]] = dol.pop(str(tmp_lst[i]))
    return dol

def peak_discovery(dol):
    #### This function works by returning the indices of possible peaks as per the definition of a local maximum
    keys = list(dol.keys())
    peaks_key = {key : [0] for key in keys}
    lop, new_pk = [], []
    ### Here the for loops over the channel for a specific nucleotide and selects possible peaks
    for key in keys:
        addendum =[50]*10
        for i in range(2, len(dol[key])):
            if dol[key][i-1] - dol[key][i-2] >= 0 and dol[key][i] - dol[key][i-1] <= 0 and dol[key][i - 1] > 0.25*sum(addendum)/10 and i not in list(range(peaks_key[key][-1], peaks_key[key][-1] + 7)): #### Wherever the derivatives switch signs there must be a local peak (min or max)
                peaks_key[key].append(i-1) ### If there is a change in derivative sign and the intensity is higher than some threshold
                addendum.append(dol[key][i-1])
                addendum.pop(0)
        peaks_key[key].pop(0)
        lop = lop + peaks_key[key]
    lop.sort()
    current = int(lop[0])
    new_d = {k : [] for k in keys}
    for j in lop[1:]: ### This is added to flag possible peaks wherever there is a substancial distance to a prior peak
        if j - current > 15 and j - current < 30:
            new_pk.append(current + 12) ### A peak position is added to the list
            pos_lst = [dol[key][current + 12] for key in keys]
            key_val = keys[pos_lst.index(max(pos_lst))]
            new_d[key_val].append(current + 12)
        elif j - current > 30: ### Same thing as before but for more spacious gaps
            new_pk.append(current + 12)
            new_pk.append(current + 26)
            pos_lst = [dol[key][current + 12] for key in keys]
            key_val = keys[pos_lst.index(max(pos_lst))] ### The peak position implied is added in the channel corresponding to the maxima
            new_d[key_val].append(current + 12)
            pos_lst = [dol[key][current + 26] for key in keys]
            key_val = keys[pos_lst.index(max(pos_lst))]
            new_d[key_val].append(current + 26)
        current = j
    for key in keys:
        peaks_key[key] = list(peaks_key[key]) + list(new_d[key])
        peaks_key[key].sort()
    lop = lop + new_pk
    lop.sort()
    seq = []
    uniq_lop_d = dict(zip(lop[:], [0]*len(lop)))
    uniq_lop = list(uniq_lop_d.keys())
    for i in uniq_lop:
        for key in list(peaks_key.keys()):
            if i in list(peaks_key[key]):
                seq.append(key)
    vals = []
    for value in peaks_key.values():
        vals = vals + value 
    return peaks_key, lop, seq #### This returns the peaks as a dictionary and as a list and the sequence of all peaks

def width(dol, dop):
    ### This function returns the width of peaks by traversing the dictionary of lists until the intensity drops under 20 units
    amplitude = {key : [] for key in list(dop.keys())}
    for key in list(dop.keys()):
        for i in dop[key]:
            k = i
            j = i
            try:
                while dol[key][i] > 20:
                    i = i + 1
                while dol[key][j] > 20:
                    j = j - 1
                amplitude[key].append(i - j)
            except IndexError:
                amplitude[key].append(max(i - k, k - j)) ### If it were to happen that the peaks have been cut and pasted an error might occur 
    return amplitude

def filterer(dol, dop):  # Here dol refers to dictionary of lists aka channels, and dop to a dictionary of peaks
    ### This function uses all previous functions to generate indentificative parameters for each peak
    length = min([len(dol[c]) for c in list(dol.keys())])
    seq_gd = list(dol.keys())
    seq_gd = "".join(seq_gd)
    dop_1 = dop
    dop_2 = copy.deepcopy(dop_1) 
    every = confidence(dop_2, dol) 
    conf = every[0]
    intens = every[1]
    peakdis_1 = every[2]
    peakdis_2 = every[3]
    amp = width(dol, dop_1)
    der_1 = {key : [] for key in list(dol.keys())}
    der_2 = {key : [] for key in list(dol.keys())}
    fullist = 0
    for key in list(dop_1.keys()): #### This must also be changed to account for channel peaks
        for i in list(dop_1[key]):
            der_1[key].append(dol[key][i] - dol[key][i - 1])
            der_2[key].append(dol[key][i + 1] - dol[key][i])
            fullist = fullist + 1
    joined, derl1, derl2, ampl = [], [], [], []
    sum_peak, sum_conf, sum_amp = [0]*10, [0]*10, [0]*10
    for i in range(0, fullist):
        try:
            vals = [dop_1[key][0] for key in list(dop_1.keys())]
        except IndexError:
            leng = [len(dop_1[key]) for key in list(dop_1.keys())]
            while len(leng) != 0 and min(leng) == 0:
                ind = leng.index(min(leng))
                der_1.pop(list(dop_1.keys())[ind], None)
                der_2.pop(list(dop_1.keys())[ind], None)
                amp.pop(list(dop_1.keys())[ind], None)
                dop_1.pop(list(dop_1.keys())[ind], None)
                leng = [len(dop_1[key]) for key in list(dop_1.keys())]
            vals = [dop_1[key][0] for key in list(dop_1.keys())]
        derl1.append(der_1[list(dop_1.keys())[vals.index(min(vals))]][0])
        derl2.append(der_2[list(dop_1.keys())[vals.index(min(vals))]][0])
        ampl.append(amp[list(dop_1.keys())[vals.index(min(vals))]][0])
        der_1[list(dop_1.keys())[vals.index(min(vals))]].pop(0)
        der_2[list(dop_1.keys())[vals.index(min(vals))]].pop(0)
        amp[list(dop_1.keys())[vals.index(min(vals))]].pop(0)
        dop_1[list(dop_1.keys())[vals.index(min(vals))]].pop(0)
        sum_peak.append(intens[i])
        sum_peak.pop(0)
        sum_conf.append(conf[i])
        sum_conf.pop(0)
        sum_amp.append(ampl[i])
        sum_amp.pop(0)
        mean_i_1, mean_c_1, mean_amp_1, mean_i_2, mean_c_2, mean_amp_2 = 0, 0, 0, 0, 0, 0
        for k in range(0,10):
            mean_amp_2 = mean_amp_2 + sum_amp[k]
            mean_c_2 = mean_c_2 + sum_conf[k]
            mean_i_2 = mean_i_2 + sum_peak[k]
            if k % 2 == 0:
                mean_i_1 = mean_i_1 + sum_peak[- int((k + 2)/2)]
                mean_c_1 = mean_c_1 + sum_conf[- int((k + 2)/2)]
                mean_amp_1 = mean_amp_1 + sum_amp[- int((k + 2)/2)]
        joined.append([conf[i], intens[i], ampl[i], mean_i_2/10, mean_amp_2/10, mean_c_2/10, mean_i_1/5, mean_c_1/5, mean_amp_1/5, peakdis_1[i], peakdis_2[i], derl1[i], derl2[i], min(vals)/length, min(vals)])
    return joined


### This can be further automated by searching for -F_ or -R_ in the path
def reader (path, fmt="abi"):
   ### This function takes the path of the forward or reverse, extracts the sequences and produces an object containing the fasta sequence
   record = SeqIO.read(path, fmt)
   c = ["DATA9", "DATA10", "DATA11", "DATA12"]
   channels = dict()
   if "-R_" in path:
      sequence = str(record.annotations["abif_raw"]["PBAS2"]).replace("b","").replace("'","")[::-1]
      sequence = sequence.replace("A","Z").replace("G","X").replace("T","A").replace("C","G").replace("Z","T").replace("X","C")
      sequence = str(">Reverse" + '\n' + sequence + '\n')
      for i in c:
         channels[i] = record.annotations["abif_raw"][i][::-1]
      # ploc = tuple(np.subtract(list(itertools.repeat(len(channels["DATA9"]), len(record.annotations["abif_raw"]["PLOC2"]))), list(record.annotations["abif_raw"]["PLOC2"][::-1])))
      ploc = [-(j - len(channels["DATA9"])) for j in record.annotations["abif_raw"]["PLOC2"][::-1]]
      guide = str(record.annotations["abif_raw"]["FWO_1"]).replace("b","").replace("'", "").replace("A","Z").replace("G","X").replace("T","A").replace("C","G").replace("Z","T").replace("X","C")

   elif "-F_" in path:
      sequence = str(record.annotations["abif_raw"]["PBAS2"]).replace("b","").replace("'","")
      sequence = str(">Forward" + '\n' + sequence +'\n')
      for i in c:
         channels[i] = record.annotations["abif_raw"][i] 
      ploc = list(record.annotations["abif_raw"]["PLOC2"])
      guide = str(record.annotations["abif_raw"]["FWO_1"]).replace("b","").replace("'", "")
   return sequence, ploc, channels, guide



def main(path_fw, path_rv):
   try:
      name=path_fw.split("\\")[-1].split("_")[2] + ".txt"
      file=open(name, "x")
   except IOError:
      file=open(name, "w")
   file.write("BASECALLER_SEQUENCES" + '\n' + str(reader(path_fw)[0]+reader(path_rv)[0]) + '\n\n\n\n\n')
   file.write('\n' + "PEAK_LOC_FW" + '\n' + str(reader(path_fw)[1]))
   file.write('\n' + "PEAK_LOC_RV" + '\n' + str(reader(path_rv)[1]))
   file.write('\n\n\n\n\n')
   file.write('\n' + "CHANNELS_FW" + '\n' + str(reader(path_fw)[2]))
   file.write('\n' + "CHANNELS_RV" + '\n' + str(reader(path_rv)[2]))
   file.write('\n\n\n\n\n')
   file.write('\n' + "GUIDE_FW" + '\n' + str(reader(path_fw)[3]))
   file.write('\n' + "GUIDE_RV" + '\n' + str(reader(path_rv)[3]))
   file.write('\n\n\n\n\n')
   file.close()
   return