import os
import subprocess
from Bio import SeqIO
import copy

#### PERHAPS A CONFIDENCE METRIC CAN BE CALCULATED AS THE VALUE OF EACH CHANNEL DIVIDED BY THE TOTAL AND CHOOSE ONWARDS FROM A THRESHOLD OF HIGH REPETITIONS
def confidence(dop, dol): 
    conl = {key : [] for key in list(dop.keys())}
    intens = {key : [] for key in list(dop.keys())}
    conf = []
    ints = []
    peaks = []
    peak_dist_fw = []
    peak_dist_bw = []
    for keys in list(dop.keys()):
        for i in dop[keys]:
            peaks.append(i)
            sums = sum([dol[c][i] for c in dol.keys()])
            intens[keys].append(dol[keys][i])
            try:
                conme = dol[keys][i]/sums
            except ZeroDivisionError:
                conme = 0
            conl[keys].append(conme)
    peaks.sort()
    for i in range(0, len(peaks)):
        try:
            peak_dist_fw.append(peaks[i + 1] - peaks [i])
        except IndexError:
            peak_dist_fw.append(-1)
        try:
            peak_dist_bw.append(peaks [i] - peaks[i - 1]) #### Value at [-1] exists
        except i == 0:
            peak_dist_bw.append(-1)
        try:
            vals = [dop[key][0] for key in list(dop.keys())]
        except IndexError:
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
    tmp_lst = list(dol.keys())
    for i in range(0, len(dol.keys())):
        dol[guide[i]] = dol.pop(str(tmp_lst[i]))
    return dol

def peak_discovery(dol):
    #### This function works by returning the indices of possible peaks as per the definition of a local maximum
    keys = list(dol.keys())
    peaks_key = {key : [0] for key in keys}
    lop = []
    for key in keys:
        for i in range(2, len(dol[key])):
            if i not in range(peaks_key[key][-1] - 3, peaks_key[key][-1] + 3) and dol[key][i-1] - dol[key][i-2] >= 0 and dol[key][i] - dol[key][i-1] <= 0 and dol[key][i - 1] > 50: #### Wherever the derivatives switch signs there must be a local peak (min or max)
                peaks_key[key].append(i-1)
        lop = lop + peaks_key[key]
    lop.sort()
    current = int(lop[0])
    new_pk = []
    new_d = {k : [] for k in keys}
    for j in lop[1:]:
        if j - current > 15 and j - current < 30:
            new_pk.append(current + 12)
            pos_lst = [dol[key][current + 12] for key in keys]
            key_val = keys[pos_lst.index(max(pos_lst))]
            new_d[key_val].append(current + 12)
        elif j - current > 30:
            new_pk.append(current + 12)
            new_pk.append(current + 26)
            pos_lst = [dol[key][current + 12] for key in keys]
            key_val = keys[pos_lst.index(max(pos_lst))]
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
    return peaks_key, lop, seq #### This returns the peaks as a dictionary and as a list

def jiggler(lop, dol): #### This function will return a list of jiggled peak locations as severance that they are locally peaks
    newd = {key : [] for key in list(dol.keys())}
    new = []
    for i in lop:
        for c in dol.keys():
            j = i - 3
            success = False
            while j < i + 3:
                if dol[c][j] >= dol[c][j + 1] and dol[c][j] > 50 and dol[c][j + 1] > 50:
                    newd[c].append(j)
                    new.append(j)
                    success = True
                    break
                elif dol[c][j] > 50 and dol[c][j + 1] > 50: 
                    j = j + 1
                    success = False
                else:
                    success = None
                    break
            if success == False:
                newd[c].append(i)
                new.append(i)
    new_sort = {key : 0 for key in new}
    newl = list(new_sort.keys())
    newl.sort()
    return newd, newl
#### This lines of code cross check the discovered peaks and the ploc from the basecaller 
def cross_check(dol, lop):
    all_peaks = peak_discovery(dol)[1]
    new_peaks = all_peaks.copy()
    ploc_fw = jiggler(lop, dol)[1]
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

def width(dol, dop): #### MUST BE CHANGED AS WELL
    amplitude = {key : [] for key in list(dop.keys())}
    for key in list(dop.keys()):
        for i in dop[key]:
            k = i
            j = i
            try:
                while dol[key][i] > 50:
                    i = i + 1
                while dol[key][j] > 50:
                    j = j - 1
                amplitude[key].append(i - j)
            except IndexError:
                amplitude[key].append(max(i - k, k - j))
    return amplitude

def filterer(dol, dop):  #Here dol refers to dictionary of lists aka channels, and dop to a dictionary of peaks
    length = min([len(dol[c]) for c in list(dol.keys())])
    seq_gd = list(dol.keys())
    seq_gd = "".join(seq_gd)
    dop_1 = dop
    dop_2 = copy.deepcopy(dop_1) #### Here the jiggler should return a dop a dictionary of peaks CHECK
    every = confidence(dop_2, dol) #### This should take into account the channels of the peaks CHECK
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
    joined = []
    derl1 = []
    derl2 = []
    ampl = []
    sum_peak = [0]*10
    sum_conf = [0]*10
    sum_amp = [0]*10
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
        mean_i_1 = 0
        mean_c_1 = 0
        mean_amp_1 = 0
        mean_i_2 = 0
        mean_c_2 = 0
        mean_amp_2 = 0
        for k in range(0,10):
            mean_amp_2 = mean_amp_2 + sum_amp[k]
            mean_c_2 = mean_c_2 + sum_conf[k]
            mean_i_2 = mean_i_2 + sum_peak[k]
            if k % 2 == 0:
                mean_i_1 = mean_i_1 + sum_peak[- int((k + 2)/2)]
                mean_c_1 = mean_c_1 + sum_conf[- int((k + 2)/2)]
                mean_amp_1 = mean_amp_1 + sum_amp[- int((k + 2)/2)]
        if seq_gd == 'CTAG':
            joined.append([conf[i], intens[i], ampl[i], mean_i_2/10, mean_amp_2/10, mean_c_2/10, mean_i_1/5, mean_c_1/5, mean_amp_1/5, peakdis_1[i], peakdis_2[i], derl1[i], derl2[i], (length - min(vals))/length, length - min(vals)]) ## Here the value for the position in the reverse must be flipped
        elif seq_gd == 'GATC':
            joined.append([conf[i], intens[i], ampl[i], mean_i_2/10, mean_amp_2/10, mean_c_2/10, mean_i_1/5, mean_c_1/5, mean_amp_1/5, peakdis_1[i], peakdis_2[i], derl1[i], derl2[i], min(vals)/length, min(vals)])
    return joined
    
    
def tmp_rm():
   for i in range(0, len(os.listdir())):
      check=str("aln_tmp_" + str(i))
      if os.path.exists(check):
         os.remove(check)
   return 

def tmp_aln(seq):
   ### This function creates a file with both sequences
   files=os.listdir()
   for i in range(0, len(os.listdir())):
      check=f"aln_tmp_{i}"
      if check in files:
         current_tmp = open(check, "r")
         line = current_tmp.read()
         if seq in line:
            current_tmp.close()
            break
         else:
            continue
      else:
         new_tmp=open(check, 'x')
         new_tmp.write(seq)
         new_tmp.close()
         break
   return check


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


def aligner(file):
   cmd = [r"C:\Mafft\mafft.bat", "--clustalout", "--localpair", "--maxiterate", "1000", "--op", "2.0", "--ep", "0.1", file]
   process = subprocess.run(cmd, capture_output=True, text=True)
   matches = process.stdout.split('\n')
   fw_seq=""
   rv_seq=""
   match_seq=""
   for i in matches:
      if "Forward" in i:
         fw_seq=fw_seq+(i[16:])
      elif "Reverse" in i:
         rv_seq=rv_seq+(i[16:])
      # A better way to get the alignment must be though of
      elif len(i)>50:
         match_seq=match_seq+(i[16:])
   match_seq=match_seq.replace(" ", "_")
   return match_seq, fw_seq, rv_seq


def locator(path_fw, path_rv):
   all=aligner(tmp_aln(reader(path_fw)[0] + reader(path_rv)[0]))
   align=all[0]
   fw=all[1]
   rv=all[2]
   Forward=reader(path_fw)[0].removeprefix(">Forward\n")
   Reverse=reader(path_rv)[0].removeprefix(">Reverse\n")
   i=0
   while i < (len(align)-20):
      kmer_fw=align[i:i+20]
      if kmer_fw.count('*') > 10:
         pos_1=Forward.index(fw[i + kmer_fw.index('*'):i + 25].replace("-","").upper())
         pos_2=Reverse.index(rv[i + kmer_fw.index('*'):i + 25].replace("-","").upper())
         break
      else:
         i=i+1
   aln_fw_ind = i
   i = -21
   while i > (-len(align)+20):
      kmer_bw=align[i:i+20]
      if kmer_bw.count('*') > 10:
         pos_3=Forward.index(fw[i - 30 : i - (len(kmer_bw) - kmer_bw.index('*'))].replace("-","").upper())
         pos_4=Reverse.index(rv[i - 30 : i - (len(kmer_bw) -  kmer_bw.index('*'))].replace("-","").upper())
         break
      else:
         i=i-1
   aln_rv_ind = len(align) + i
   return [pos_1, pos_3], [pos_2, pos_4], aln_fw_ind, aln_rv_ind


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
   file.write('\n' + "ALIGNMENT" + '\n' + aligner(tmp_aln(reader(path_fw)[0] + reader(path_rv)[0]))[0])
   file.write('\n\n\n\n\n')
   file.write('\n' + "ALIGNED_FW" + '\n' + aligner(tmp_aln(reader(path_fw)[0] + reader(path_rv)[0]))[1])
   file.write('\n\n\n\n\n')
   file.write('\n' + "ALIGNED_RV" + '\n' + aligner(tmp_aln(reader(path_fw)[0] + reader(path_rv)[0]))[2])
   file.write('\n\n\n\n\n')
   file.write('\n' + "LOCATORS" + '\n' + str(locator(path_fw, path_rv)))
   file.close()
   return