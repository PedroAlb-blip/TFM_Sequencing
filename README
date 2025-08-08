PREREQUISITES 

INTRODUCTION

Hi fellow user! This readme file will hopefully serve you as both a manual and an aid in the installation and proper use of the python code files present
in the master branch. All these files are a part of my master's thesis and they are supposed to produce a consensus sequence from forward and reverse files
encoded in .ab1 files produced by Sanger sequencing or any other kind of sequencing. This file will be written in a more narrative and entertaining manner
so as to improve the ease of reading and to prove that bioinformaticians are NOT boring. 

To start from the top, whilst sequencing two types of file are produced: .seq files and .ab1 files. The former contain the sequence in raw text as produced
by the sequencer written always from 5' to 3', problem is that these files tend to be quite unreadable because the base caller, i.e. the algorithm within
the sequencer used to produce the sequence has less initiative than an introvert on a friday night. What I mean is that whenever some nucleotide might have
ambiguity within the electropherogram, i.e. the chromatogram of fluorescence detection for the nucleotides produced by sequencing, the base caller will
completely forgo base calling and insert Ns, ambiguous nucleotides, in their stead making them of not much interest. The latter files, on the other hand,
are like the Swiss knife of sequencing, basically every single detail of the sequencing process is projected within these files: electrohperogram data in
various stages of denoising and normalisation, the .seq file sequence, the peak location or ploc of each of the nucleotides, even when ambiguous, within
the electropherogram data, values for confidence to produce fastq sequences wherever they may apply, even metadata regarding the start and end times of the
sequencing process, the name of the person who sent the order for sequencing, their social security number, mode of payment and latest internet browser
history. 

Well perhaps not the latter three, but you get the gist, everything and anything needed to reproduce the electropherograms is witin the .ab1 files (the
everything files), however handling these files has a big drawback, they are not produced as raw text but rather as binary code which is only readable to a
number of programs like BioEdit or Mega, luckily multiple packages in various programming languages exist which can allow us to read and extract data from
the .ab1 files to then redo the peak detection, and try to generally speaking produce not so shitty sequences which can be subject to alignent. And thus it
begins, the endless struggle to produce higher and higher quality of reads, faster and more precisely than could otherwise be achieved manually. 

PARSING

The first file generated throughout this project and, coincidentally, the first one to have to be executed is the parser.py file, the idea here is simple,
use the SeqIO package to read the .ab1 files, pick and choose relevant information and produce parsed files. Again, from the top, multiple user inputs are
needed within this file, perhaps not the most clever way of doing it but still the simplest. This file needs to be modified to add:
 > The full file path containing the .ab1 files which you want parsed the syntax here is important the home directory C has to be lower case (so c) and the
directories (folders) have to be spaced by 2 '/', e.g. "c://Home//User//Downloads//Sequence_Folder".
 > The forward and reverse primer abbreviations written within the .ab1 files, these are important to properly classify the files as forward and reverse,
as they will be treated differently as the reverse file must be reverse complemented because it refers to the sequence of the complementary strand.
 > The file names which the parsed files will be given, by default the files produced are given the name within the .ab1 file which should, emphasis here,
be the name of the sample given the moment the sequencing task is generated. 

I assume that file naming is not homogenuous everywhere and therefore this file may need some tweaking to achieve desired results. This subprogram
essentiales picks and chooses the channels, i.e. the values for fluorescence detection across time in the different fluorescence colors which have been
normalised and denoised, corresponding to DATA9-DATA12 as well as the guides, that is, a string containing the nucleotide to channel relation for both
forward and reverse. Sequences in raw text are also inserted for forward and reverse, as well as the PLOC or peak locations, the indices of peak positions
within the channels, even though these are not needed for future processing. A function called rename has also been introduced to change the names of the
channels to the corresponding nucleotide channel.

PEAK DISCOVERY

After parsing the peaks are discovered, this is performed by checking the derivatives in every point of the electropherogram at each of the channels,
exceptuating peaks near each other and intensity values under 50 intensity units. Also whenever a peak has not been detected for over 15 positions a new 
tentative peak in flagged at the channel at maximum intensity at said point. Peak discovery is later used within the filterer function so as to generate 
several parameters of perhaps key value in correctly ascribing peak status, these include confidence, i.e. the value of the peaking channel over the value 
of the total sum of intensities across all channels, the amplitude of peak, and peak intensity, as well as means for each of them with the peaks at 5 and 
10 positions backwards, derivative values pre-peak and after-peak, and distance to the next and prior peak, as well as the absolute and relative distance 
from the 5' start. 

The idea here is to train an ML model to learn from peak description introduced through manual entry, assessing peak values in a binary fashion thorugh 
graphing and manual input. 

ACME-ING

Now ACME.py, standing not for the American Company that Makes Everything from coyote and the roadrunner but rather Acquisition of Crests by Manual Entry, 
is capable of reading the parsed files, producing a graphical representation of 50 tentative peaks at a time allowing the user to input values either 0 if 
the peak is ascribed correctly, and 1 if it is not, shortcuts for all peaks considered right and none considered right were produced by typing ALL and NONE 
respectively. 

Once the binary vectors are produced they are stored within the results.txt and the parameters produced are stored within the training_params.txt these are
used within the train_valid.py, which implements a 3 dense layer keras neural network the parameters are used as input of a 15x32 input layer which is 
activated by relu, and passed onto a second hidden layer 32x16 with relu activation, finally a 16x1 output layer with sigmoid activation is utilised. In 
order to optimise weights the binary focal crossentropy loss function is used and the adamax optimiser is the one utilised. The reason behind this is that 
the output of the neural network is a binary and ideally binary focal crossentropy takes less into account easy to classify peaks while further penalising 
misclassification. And adamax is also a better optimiser in this instance, this was found by trial and error but corroborated by a quick literature 
reasearch as it is less susceptible to vanishing gradients than the ADAM optimiser.

The entire model is saved within a .keras file called checkpoint.keras and the better model is saved as checkpoint.keras, and can be loaded later by a 
keras.model.load(). Also loss is graphed across epochs, as well as accuracy, these two parameters are used to determine whether the model has properly 
learnt or not, if you are to use the model introduced, beware only 95% accuracy was achieved, however after some sequence processing it would seem that 
only very noisy sequences are affected by the subpar accuracy, which will allow me to sleep soundly at night. 

DiNNA-ING

The DiNNA.py file is capable of taking the 
