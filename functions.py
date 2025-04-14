import os
import subprocess
from Bio import SeqIO

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
                if dol[c][j] >= dol[c][j+1] and dol[c][j] > 50 and dol[c][j + 1] > 50:
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