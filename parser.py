##### This programs are supposed to parse the ab1 forward and reverse files produced by sanger sequencing 
##### The output of these programs is a consensus sequence between forward and reverse which will be attained using the following steps
##### 1.- The ab1 files will be parsed to extract sequence information (PBAS2), peak location (PLOC2), and fluorescence values and channels (FWO_1, DATA9 through 12)
##### 2.- The forward and reverse sequences will be transplanted onto a temporary file and an alignment be made
##### 3.- From said alignment the position of the first match will be mapped and stored 
##### 4.- Perhaps a new file containing sequences, PLOC, alignment and match position can be produced to save time and effort for later usage

import os
from Bio import SeqIO


### This can be further automated by searching for -F_ or -R_ in the path this varies between file names
def reader (path, fmt="abi"):
   ### This function takes the path of the forward or reverse, extracts the sequences and produces an object containing the fasta sequence
   ### Therefore the indicator for fw and rv must be introduced 
   record = SeqIO.read(path, fmt)
   c = ["DATA9", "DATA10", "DATA11", "DATA12"]
   channels = dict()
   if "-R_" in path: ### Here the indicator for reverse primer
      sequence = str(record.annotations["abif_raw"]["PBAS2"]).replace("b","").replace("'","")[::-1]
      sequence = sequence.replace("A","Z").replace("G","X").replace("T","A").replace("C","G").replace("Z","T").replace("X","C")
      sequence = str(">Reverse" + '\n' + sequence + '\n')
      for i in c:
         channels[i] = record.annotations["abif_raw"][i][::-1]
      # ploc = tuple(np.subtract(list(itertools.repeat(len(channels["DATA9"]), len(record.annotations["abif_raw"]["PLOC2"]))), list(record.annotations["abif_raw"]["PLOC2"][::-1])))
      ploc = [-(j - len(channels["DATA9"])) for j in record.annotations["abif_raw"]["PLOC2"][::-1]]
      guide = str(record.annotations["abif_raw"]["FWO_1"]).replace("b","").replace("'", "").replace("A","Z").replace("G","X").replace("T","A").replace("C","G").replace("Z","T").replace("X","C")

   elif "-F_" in path: ### Here the indicator for forward primer
      sequence = str(record.annotations["abif_raw"]["PBAS2"]).replace("b","").replace("'","")
      sequence = str(">Forward" + '\n' + sequence +'\n')
      for i in c:
         channels[i] = record.annotations["abif_raw"][i] 
      ploc = list(record.annotations["abif_raw"]["PLOC2"])
      guide = str(record.annotations["abif_raw"]["FWO_1"]).replace("b","").replace("'", "")
   return sequence, ploc, channels, guide



def main(path_fw, path_rv):
   ### This function generates a new file containing the name of the sample gotten by processng of the string of the file name
   ### As well as the channels, plocs, sequences, and guides for each of them
   try:
      name=path_fw.split("\\")[-1].split("_")[2] + ".txt"
      file=open(f"Sequences\\{name}", "x")
   except IOError:
      file=open(f"Sequences\\{name}", "w")
   file.write("BASECALLER_SEQUENCES" + '\n' + str(reader(path_fw)[0]+reader(path_rv)[0]) + '\n\n\n\n\n')
   file.write('\n' + "PEAK_LOC_FW" + '\n' + str(reader(path_fw)[1]))
   file.write('\n' + "PEAK_LOC_RV" + '\n' + str(reader(path_rv)[1]))
   file.write('\n\n\n\n\n')
   file.write('\n' + "CHANNELS_FW" + '\n' + str(reader(path_fw)[2]))
   file.write('\n' + "CHANNELS_RV" + '\n' + str(reader(path_rv)[2]))
   file.write('\n\n\n\n\n')
   file.write('\n' + "GUIDE_FW" + '\n' + str(reader(path_fw)[3]))
   file.write('\n' + "GUIDE_RV" + '\n' + str(reader(path_rv)[3]))
   file.close()
   return

##### VERY IMPORTANT: POSITIONS OF PLOC CHANGE RELATIVE TO THE SEQUENCE WHEN DOING A REVERSE COMPLEMENT 
##### Perchance doing the absolute value of the subtraction between the PLOCs and the total length 14000


path_folder = "c:\\Users\\Pedro\\Downloads\\resecuenciasdesp10asp29" ### Here the full path to the folder containing the .ab1 files
file_list = [j for j in os.listdir(path_folder) if ".ab1" in j]
files = [file.split("_")[2] for file in file_list]
file_purge = dict(enumerate(files))
egrup_elif = {v : k for k, v in file_purge.items()}
unique = list(egrup_elif.keys())
for i in unique:
   fw, rv = "", ""
   for j in file_list:
      if i in j and "-F_" in j: ### Here the indicator for forward primer
         fw = str(path_folder + "\\" + j)
      elif i in j and "-R_" in j: ### Here the indicator for reverse primer
         rv = str(path_folder + "\\" + j)
   if rv != "" and fw != "":
      main(fw, rv)
