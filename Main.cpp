/*
Developed by: Michele Busby, Computational Biologist
Broad Technology Labs/MBRD
Broad Institute 
Last Update: 4/27/16

/*=================================
 Useage: 
  ./GetRandomByReadNormalizedByInsertSize -bam bam1.bam -bam bam2.bam -bam bam3.bam -out /outDirectory/ 
 
 
=================================*/


#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <set>

#include "api/BamReader.h"
#include "api/BamWriter.h"
#include "Handy.h"

//g++ *.cpp -L/seq/mbrd/mbusby/Software/bamtools/lib -I/seq/mbrd/mbusby/Software/bamtools/include -I/seq/mbrd/mbusby/Software/bamtools/include/api



unsigned int checkErrors();
void displayHelp();
void somethingsGoneWrong(string);
void processBamFiles();
void getFragmentSizeMaps();
inline char separator();

//Variables taken from the inputs

string outDir="";
int seed=32413;
bool proportional=false;

vector <string> bamFileNames;
string controlBam="";
bool hasControl=false;

map<string, double> randMap;

vector < map<int, int> > fragmentSizes;

vector < map<int, double> > samplingProbabilities; 

void checkBamsForOpening();
void getFragmentSizeMaps();
void getProbabilities();
void getProbabilitiesForControl(map<int, int>);
int getTotalFromIntIntMap(map<int,int>);
bool is_folder_writable(const char* );


using namespace std;
using namespace BamTools;


int main(int argc, char* argv[]) 
{
		
	//Intialize library of useful things
	Handy h(0);
	
	int optind=1;

	while ((optind < argc) && (argv[optind][0]=='-')) 
	{	
        string sw = argv[optind];
				
		if (sw=="-h") 
		{	
            optind++;
			displayHelp();
			return 1;
        }
		
		else if(optind >=  argc-1)
		{
			cerr<<"Your final parameter, "<<sw<<" is missing a value."<<endl;
			return 1;
		}

		else if (h.cmpStringNoCase(sw, "-bam")==1)
		{	
            optind++;
			bamFileNames.push_back(argv[optind]);			
			optind++;
        }
		
		else if (h.cmpStringNoCase(sw, "-control_bam")==1)
		{	
            optind++;
			controlBam=(argv[optind]);			
			optind++;
        }
		
		else if (h.cmpStringNoCase(sw, "-out_dir")==1)
		{	
            optind++;
			outDir = argv[optind];		
			optind++;
        }	
		
		
		else if (h.cmpStringNoCase(sw, "-seed")==1)
		{	
            optind++;
			seed=h.getIntFromString(argv[optind]);
			optind++;
		
        }
		
		else if (h.cmpStringNoCase(sw, "-mode")==1)
		{	
            optind++;
			if(h.cmpStringNoCase(argv[optind], "proportional"))
			{
				proportional=true;
			}
			
			optind++;
		
        }
		
		
		else
		{
			cerr<<"Main: Unknown parameter:"<<sw<<endl;
			return 1;
		}
	}	
	
	cout<<"Beginning program. Assuming reads are uniquely aligned and deduplicated!"<<endl;
	cout<<"Assuming all bam files have different names or they will overwrited each other!"<<endl;
	cout<<"Assuming all reads are the same length (first and second reads can differ))"<<endl;
	
	cout<<"Bam files to be normalized:"<<endl;
	
	for(int i=0; i<bamFileNames.size(); ++i)
	{
		cout<<"    "<<bamFileNames[i]<<endl;		
	}
	
	if(controlBam.length()>0)
	{		
		cout<<"Control bam:"<<controlBam<<endl;	
		hasControl=true;
		bamFileNames.push_back(controlBam);
	}
	
	checkBamsForOpening();
	
	checkErrors();

	getFragmentSizeMaps();
	
	getProbabilities();
	
	processBamFiles();
	
	cout<<"Done\n";
		
}

/*************************************************
Outputs the sampled bam

*************************************************/
void processBamFiles()
{
	Handy h(0);
	BamReader reader;
	BamAlignment al;
	RefVector refVector;
	double X;

	srand (seed);	
	
	string outFileName;
	
	
	
	for(unsigned int i=0; i<bamFileNames.size(); ++i)
	{
		string bamFileName=bamFileNames[i];
		
		cout<<"Processing "<<bamFileName<<endl;
		
		cout<<"Size of sampling probabilities:"<<samplingProbabilities.size()<<endl;
		map<int, double> probabilities=samplingProbabilities[i];
		
		
		set<string> readsToWrite;
		
		//Open writer (kind of early in case things fail but not so early as to not be super annoying)
		std::string base_filename = bamFileName.substr(bamFileName.find_last_of(separator()) + 1);		
		
		outFileName=outDir+separator()+base_filename+".normalized.bam";
		
		
		cout<<outFileName<<endl;		
				
		ofstream outputStream;
		outputStream.open(outFileName.c_str());		
		
		
		//Second iteration through bam - 
		//Coin flip to make a map of all the reads the program is going to keep
	
		if ( !reader.Open(bamFileName) ) {
			cerr << "Could not open input BAM file" << bamFileName << endl;
			return;
		}
					
		
		cout<<"Performing second of three iterations through "<<bamFileName<<" Coin flip. "<<endl;
	
		while(reader.GetNextAlignment(al))
		{		
			
			//Only operates on first mate so that both pairs will be in there
			//Assumes both mates are in bam with insert size =0 when one is not mapped
			if( al.IsFirstMate()==true )
			{
				X=((double)rand()/(double)RAND_MAX);
				
				//where the probaility table has the probability of keeping the read
				if(X<=probabilities[abs(al.InsertSize)] )
				{
					readsToWrite.insert(al.Name);
				}
			}
		}			
		
		reader.Close();
		
		//Third iteration through bam 
		//Closes and starts over again at the beginning. There's probably a way to rewind rather than doing this but this works.
		//Check that the read won the coin flip and writes the bams it is going to keep
		
		cout<<"Performing third of three iterations through "<<bamFileName<<" Writing reads that won the coin flip to output bam. "<<endl;
		
		if ( !reader.Open(bamFileName) ) {
			cerr << "Could not open input BAM file the second time, which is weird: " << bamFileName << endl;
			return;
		}
		
		//Set up Writer and grab the headers and paste them in the new output files.
		BamWriter::CompressionMode compressionMode = BamWriter::Compressed;
		
		// open BamWriter
		BamWriter writer;	
		
		// get BamReader metadata
		const string headerText = reader.GetHeaderText();
		const RefVector references = reader.GetReferenceData();
		if ( references.empty() ) {
			cout << "bamtools random ERROR: no reference data available...Trying to carry on." << endl;        
		}
		
		writer.SetCompressionMode(compressionMode);
		if ( !writer.Open(outFileName, headerText, references) ) {
			cerr << "bamtools random ERROR: could not open " << outFileName
				 << " for writing... Aborting." << endl;
			reader.Close();
			return;
		}
		
		cout<<"Writing alignments to "<<outFileName<<endl;
		
		//Write both pairs for all the reads that won the coin flip during the second iteration
		while(reader.GetNextAlignment(al))
		{				
			
			if(readsToWrite.find(al.Name) != readsToWrite.end() )
			{
				writer.SaveAlignment(al);				
			}
		}	
	
		
		reader.Close();
		writer.Close();	
		//Make sure the set gets cleaned up before it uses to much memory and crashes the thing
		readsToWrite.clear();
		
	}
		
}

/*=============================================================
Get map of fragment sizes where fragment size = firstReadLength + insert size + second read length
That should give you the size of the DNA fragment that got sequenced
===============================================================*/
void getFragmentSizeMaps()
{
	Handy h(0);
	BamReader reader;
	BamAlignment al;
	RefVector refVector;
	string readName;
	
	for (unsigned i=0; i < bamFileNames.size(); i++) {	

		map <int , int  > fragSz;
		
		string bamFileName=bamFileNames[i];
		
		
		
		if ( !reader.Open(bamFileName) ) 
		{
			cerr << "Could not open input BAM file." << endl;
			return;
		}
		
		unsigned int ctr=0;
			
		//Count the distribution of fragment sizes for the bams
		while(reader.GetNextAlignment(al))
		{	
			++ctr;		
			if(ctr%10000000==0)
			{
				cout<<"Reading line "<<ctr<<" of file "<<i<<" "<< bamFileName <<" for insert size map."<<endl;
			}
				
			if(al.IsFirstMate() )
			{
				++fragSz[ abs(al.InsertSize) ];				
			}	
		}
		
		fragmentSizes.push_back(fragSz);
		reader.Close();
	}
}	
	

void getProbabilities()
{
	map<int, int> lowestReadCountByinsertSize;
	
	int nRegularBams=fragmentSizes.size();
	if(hasControl)
	{
		nRegularBams=nRegularBams-1;
	}
	
	
	//Get the lowest observed read count by insert size in the regular bams
	//This only passes over the regular bams so it will not be limited by the size of the WCE
	for(unsigned int i=0; i<nRegularBams; ++i)
	{
		typedef std::map<int, int>::iterator it_type;
		
		for(it_type iterator = fragmentSizes[i].begin(); iterator != fragmentSizes[i].end(); iterator++) 
		{
			
			if ( lowestReadCountByinsertSize.find(iterator->first) == lowestReadCountByinsertSize.end() ) 
			{
				lowestReadCountByinsertSize[iterator->first]=iterator->second;
			} 
			else 
			{
				if(lowestReadCountByinsertSize[iterator->first]>iterator->second)
				{
					lowestReadCountByinsertSize[iterator->first]=iterator->second;
				}
			}
			
		}
	}
	
	//If the insert size is unobserved in any of the other bam files OR the control set it to zero
	typedef std::map<int, int>::iterator it_type;
	for(it_type iterator = lowestReadCountByinsertSize.begin(); iterator != lowestReadCountByinsertSize.end(); iterator++) 
	{		
		for(unsigned int i=0; i<fragmentSizes.size(); ++i)
		{			
			if ( fragmentSizes[i].find(iterator->first) == fragmentSizes[i].end() ) 
			{			
				lowestReadCountByinsertSize[iterator->first]=0;
			}				
		}		
	}	
	
	//Now you have the distribution for each, and the distribution of the lowest count in each bam with zeros for places where they do not appear in any one of the bams.
	//includes zeros as these are noise, bad mappings etc. Leave them in as they will represent the noise reads in each sample, i.e. some lowest X number of unpaired
	
	
	//Interate through insert sizes
	//Then iterate through map of new insert sizes
	cout<<"Number of regular bams"<<nRegularBams<<endl;
	
	for(unsigned int i=0; i<nRegularBams; ++i)
	{	
		map <int, double> probabilities;
		//For each insert size create a new map of the probabilities
		for(it_type iterator = fragmentSizes[i].begin(); iterator != fragmentSizes[i].end(); iterator++) 
		{		
			//NANs are in here but won't matter because we won't encounter that value back in the bam file
			probabilities[iterator->first]=(double)lowestReadCountByinsertSize[iterator->first]/(double)iterator->second;			
		}	
	
		samplingProbabilities.push_back(probabilities);		
	}	
		
	cout<<"All done building sampling probabilities for regular bams!"<<endl;
	cout<<"Beginning sampling probability for control!"<<endl;
	
	if(hasControl)
	{
		getProbabilitiesForControl(lowestReadCountByinsertSize);
	}
}	
	
	
/*=========================================================================================
Get proportional probability tables for Control
Got: fragmentsizes,
 lowest read count by insert size for each of them except the control
 
Logic:
	You can't add reads so you want to cut off as few as possible
	You want to cut off reads in the control rather than the bams because the control contains less information
	There will be one fragment size that limits the distribution	
	For each fragment size, get the total number of possible reads if that is the proportion of the reads
	What is the lowest number of possible reads for the control bam for each fragment size (how much do you have to cut off to make it fit)
 
===========================================================================================*/
void getProbabilitiesForControl(map<int, int> lowestReadCountByinsertSize){
	
	map<int, int> controlCt=fragmentSizes.back();
	map<int, double> pdfLowCounts;
	
	int totalReads=getTotalFromIntIntMap(lowestReadCountByinsertSize);
	int totalControlReads=getTotalFromIntIntMap(controlCt);
	cout<<"totalControlReads"<<totalControlReads<<endl;
	
	//Get the proportion of each for the lowest count distribution
	typedef std::map<int, int>::iterator it_type;
	for(it_type iterator = lowestReadCountByinsertSize.begin(); iterator != lowestReadCountByinsertSize.end(); iterator++) 
	{
		pdfLowCounts[iterator->first]=(double) iterator->second/(double)totalReads;
		
	}
	
	//total reads something less than or equal to total control reads
	int minTotalReads=totalControlReads;
	cout<<"Min total reads:"<<minTotalReads<<endl;
	
	for(it_type iterator = controlCt.begin(); iterator != controlCt.end(); iterator++) 
	{	
		
		//read count / portion of total reads = total reads needed for this to be true
		int thisMin=floor( (double) iterator->second /(double)pdfLowCounts[iterator->first]);
		
		if(thisMin>0)
		{
			minTotalReads=min(minTotalReads,thisMin);
		}
		cout<<iterator->first<<" "<<iterator->second<<" "<<minTotalReads<<endl;
	}
	
	
	map <int, double> probabilities;
	
	typedef std::map<int, double>::iterator it_double_type;
	//For each insert size create a new map of the probabilities
	for(it_double_type iterator = pdfLowCounts.begin(); iterator != pdfLowCounts.end(); iterator++) 	
	{		
		int totalReadsThisinsertSize=floor( pdfLowCounts[iterator->first] * (double) minTotalReads);
		
		probabilities[iterator->first]=(double)totalReadsThisinsertSize/(double) controlCt[iterator->first];	//And this gives the probability you want to pull it in the coin toss	
	}	
	
	samplingProbabilities.push_back(probabilities);		
		
	
}

/*===========================================================================================
GetTheTotalReads in the map
===========================================================================================*/

int getTotalFromIntIntMap(map<int,int> intIntMap){
	
	int total=0;
	
	typedef std::map<int, int>::iterator it_type;
	for(it_type iterator = intIntMap.begin(); iterator != intIntMap.end(); iterator++) 
	{
		total=total+iterator->second;
	}
	
	return total;
}


/*===========================================================================================
Check that all of the necessary fields exist.
===========================================================================================*/

unsigned int checkErrors()
{
	//Errors 

	int err=0;
	Handy h(0);
	string problems="";
	
	if(bamFileNames.size()==0)
	{
		problems.append("A bam file containing the reads is needed (-bam).\n");
		++err;
	}
	
	if(outDir.length()==0)
	{		
		cout<<"OutDir after inner loop:"<<outDir<<endl;
		++err;
	}	
	
	if( is_folder_writable( outDir.c_str() ) )
	{
		cout<<"Folder is writable."	<<endl;
	}
	else
	{
		problems.append("Output directory cannot be written to.\n");
		++err;
	}
	
	/*============================================================
	Check that all the relevant files can be read/written to
	==============================================================*/
	for(unsigned int i=0; i<bamFileNames.size(); ++i)
	{
		err=err+h.checkRead(bamFileNames[i]);			
	}
	
	if(err>0)
	{
		somethingsGoneWrong(problems);
	}
	
	return err;
	
}

inline char separator()
{
#ifdef _WIN32
    return '\\';
#else
    return '/';
#endif
}

/*===========================================================================================
Checks that all the bams can be opened
===========================================================================================*/

void checkBamsForOpening()
{
	Handy h(0);

	
	for(unsigned int i=0; i<bamFileNames.size(); ++i)
	{
		BamReader reader;
		BamAlignment al;
		
		string bamFileName=bamFileNames[i];	
		cout<<bamFileName<<endl;
		
		bool firstRead=false;
		bool secondRead=false;
		
		if ( !reader.Open(bamFileName) ) {
			cerr << "Could not open input BAM file while assessing read lengths:" << bamFileName << endl;
			return;
		}
		
		
		while(reader.GetNextAlignment(al) && (firstRead==false || secondRead==false) )
		{	
	
			if( al.IsFirstMate()==true )
			{
				firstRead=true;
			}
			else
			{
				secondRead=true;
			}
		}
		
		if(firstRead==false || secondRead==false) 
		{
			cerr<<"ERROR: THE bam file "<<bamFileName<<" does not contain paried end reads"<<endl;
			return;
		}	
	
		reader.Close();
		
	}
	
	cout<<"All bam files successfully opened."<<endl;
}

//Checks that folder is writable
bool is_folder_writable(const char* str) {
    if(access(str, W_OK) == 0) {
        return true;
    } else {
        return false;   
    }
}



void displayHelp()
{

	cout<<"Usage:"<<endl;
	cout<<"./GetRandomByReadNormalizedByInsertSize -bam bam1.bam -bam bam2.bam -bam bam3.bam -out /outDirectory/ "<<endl;
	cout<<"  "<<endl;
	cout<<"Important parameters:\n";	
	cout<<"-bam Name of the bam file. The final name should be unique in the group."<<endl;
	cout<<"-control_bam A bam file downsampled to the same proportions as the other bams but is not downsampled."<<endl;
	cout<<"-out_dir The directory where all the output bams will be written."<<endl;
	//cout<<"-mode Mode is proportional or downsample (default). If proportional mode is chosen then output will have the same distribution of fragment sizes but not read count."<<endl;

	cout<<"Optionall:\n";
	cout<<"-seed A random number seed."<<endl;
}

void somethingsGoneWrong(string whatsGoneWrong)
{

	cout<<"ERROR: Something has gone horribly wrong.\n";	
	cout<<whatsGoneWrong;
	cout<<"\n";	
	cerr<<"\nPlease try again.";
	
}


