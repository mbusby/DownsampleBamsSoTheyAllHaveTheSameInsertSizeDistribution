Authors: Michele Busby (code and specification)
Broad Technology Labs

What does it do?
Takes in a group of bam files and downsamples them so that they all have ~ the same number of reads and ~ the same distribution of insert sizes.

The program first counts how many read pairs are present for each insert size for each of a set of aligned files. It then selects the lowest read count for each insert size from among the set of alignments. 

For example, if four alignments for a given antibody have one, two, three, and four million reads with an insert size of 100, all four alignments would be randomly downsampled so that the four normalized alignments each have ~one millon reads with an insert size of 100. This is performed for each insert size present in all of the alignments in the group to yield final bam files with ~ the same number of reads and insert size distribution with some minimal noise introduced to the total counts by the random sampling. This approach therefore allows for identical insert size distributions while maximizing the number of reads included in the output files. 

Why would I want to do that?
You might want to do this for ChIP-Seq as a normalization approach to compare two samples which have very different insert size distributions.

How does it work? 
Type:
 ./GetRandomByReadNormalizedByInsertSize -bam bam1.bam -bam bam2.bam -bam bam3.bam -out /outDirectory/ 

Required parameters:
  -bam Name of the bam file. The final name should be unique in the group. Enter a -bam for each bam you want to include
	-out_dir The directory where all the output bams will be written.


How to get it to compile:

Compiling the file may or may not be a big hassle because you need bamtools. Try using the binary.

First, though, copy all the files in here to a single directory on your unix server. The server should have g++ is installed (it probably is).


Bamtools

To get this to compile you will need to download Derek Barnett's bamtools API and install it in a folder somewhere. It is here: https://github.com/pezmaster31/bamtools You need the API, not the command line program though is quite useful to have. Then, after you install it, 

--> you need to edit the Makefile in this folder so that everywhere it says "FolderWhereBamToolsIs" you put in the folder where bamtools is located.

Compiling

Go to the folder where everything is installed in type "make". Ignore the warnings about Handy comparing integers. If nothing says ERROR and it makes an executable called ComplexityByStartPos you should be all set.
