##### This programs are supposed to parse the ab1 forward and reverse files produced by sanger sequencing 
##### The output of these programs is a consensus sequence between forward and reverse which will be attained using the following steps
##### 1.- The ab1 files will be parsed to extract sequence information (PBAS2), peak location (PLOC2), and fluorescence values and channels (FWO_1, DATA9 through 12)
##### 2.- The forward and reverse sequences will be transplanted onto a temporary file and an alignment be made
##### 3.- From said alignment the position of the first match will be mapped and stored 
##### 4.- Perhaps a new file containing sequences, PLOC, alignment and match position can be produced to save time and effort for later usage

import subprocess
import os
from Bio import SeqIO

def tmp_rm():
   for i in range(0, len(os.listdir())):
      check=str("aln_tmp_" + str(i))
      if os.path.exists(check):
         os.remove(check)
   return 

tmp_rm()

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
   # file.write('\n\n\n\n\n')
   # file.write('\n' + "ALIGNMENT" + '\n' + aligner(tmp_aln(reader(path_fw)[0] + reader(path_rv)[0]))[0])
   # file.write('\n\n\n\n\n')
   # file.write('\n' + "ALIGNED_FW" + '\n' + aligner(tmp_aln(reader(path_fw)[0] + reader(path_rv)[0]))[1])
   # file.write('\n\n\n\n\n')
   # file.write('\n' + "ALIGNED_RV" + '\n' + aligner(tmp_aln(reader(path_fw)[0] + reader(path_rv)[0]))[2])
   # file.write('\n\n\n\n\n')
   # file.write('\n' + "LOCATORS" + '\n' + str(locator(path_fw, path_rv)))
   file.close()
   return

##### VERY IMPORTANT: POSITIONS OF PLOC CHANGE RELATIVE TO THE SEQUENCE WHEN DOING A REVERSE COMPLEMENT 
##### Perchance doing the absolute value of the subtraction between the PLOCs and the total length 14000


path_folder = "c:\\Users\\Pedro\\Downloads\\secuenciasvp4primeralternativos"
file_list = [j for j in os.listdir(path_folder) if ".ab1" in j]
files = [file.split("_")[2] for file in file_list]
file_purge = dict(enumerate(files))
egrup_elif = {v : k for k, v in file_purge.items()}
unique = list(egrup_elif.keys())
for i in unique:
    for j in file_list:
        if i in j and "-F_" in j:
            fw = str(path_folder + "\\" + j)
        elif i in j and "-R_" in j:
            rv = str(path_folder + "\\" + j)
    main(fw, rv)
tmp_rm()