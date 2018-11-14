
#include <iostream>
#include <cstdlib>
#include <map>

#include <seqan/sequence.h>
#include <seqan/stream.h>
#include <seqan/store.h>
#include <seqan/basic.h>
#include <seqan/seq_io.h>
#include <seqan/arg_parse.h>


using namespace seqan;


typedef std::map< DnaString, int>  TKmerCounter;




StringSet<DnaString> loadGenome(auto const & filename)
{	
	// Parsing file
	CharString seqFileName = filename;
    StringSet<CharString> ids;
    StringSet<DnaString,Owner<> > seqs;
    SeqFileIn seqFileIn(toCString(seqFileName));
    readRecords(ids, seqs, seqFileIn);

    return(seqs);
  }

bool oldHaveLowComplexity(DnaString sequence, int threshold){
    unsigned l = length(sequence);
    std::unordered_map<char,int> counter;
    for(auto c: sequence){
        if(counter.count(c) !=0 ){
            counter[c]+=1;    
        }
        else{
            counter[c]=1;
        }
    }
    for(auto v:counter){
        if( v.second*100/l >= threshold){
            return(true);
        }
    }
    return(false);
}

bool haveLowComplexity(DnaString sequence, int threshold){
    // New version, using DUST2 method
    // scanning 2-mers, squaring count, discard if over limit
    
    unsigned l = length(sequence);
    std::unordered_map<std::string,int> counter;
    std::string seq = toCString(CharString(sequence));
    // reading using sliding window of 2
    for(int i =0; i < l-1; i++){
        std::string c = seq.substr(i,2);
        if(counter.count(c) !=0 ){
            counter[c]+=1;    
        }
        else{
            counter[c]=1;
        }
    }

    for(auto v:counter){
        if( v.second * v.second > threshold){
            return(true);
        }
    }
    return(false);
}

void countKmer(std::string fileName, int k, int threshold){

	StringSet<DnaString> sequences = loadGenome(fileName);
	// storage
	TKmerCounter kmerSet;

    for(int j = 0; j< length(sequences); j++)
    {
        DnaString seq = sequences[j];
    	for( int i = 0; i < length(seq) - k + 1; i++){

    		Infix<DnaString>::Type inf = infix(seq, i, i+k);
            DnaString km = DnaString(inf);

            if(not haveLowComplexity(km,threshold)){
                if ( kmerSet.count(km) !=0 ){
        			kmerSet[km] += 1;
                }
        		
        		else{
        			kmerSet[km] = 1;
        		}
            }

    	}
    }
	std::vector<std::pair< DnaString, int > > kmerVec(std::make_move_iterator(kmerSet.begin()), std::make_move_iterator(kmerSet.end()));
    std::sort(kmerVec.begin(), kmerVec.end(), [](auto x, auto y){ return x.second > y.second;} );


	for (auto elem : kmerVec)
    {
		std::cout << ">" << elem.first << "_" << elem.second  << std::endl;
        std::cout << elem.first << std::endl;
	
    }
}


int main(int argc, char const ** argv)
{
	 // Setup ArgumentParser.
    seqan::ArgumentParser parser("findAdapt");

    addArgument(parser, seqan::ArgParseArgument(
        seqan::ArgParseArgument::STRING, "FileName"));

    addOption(parser, seqan::ArgParseOption(
        "k", "k", "Size of the kmer to generate",
        seqan::ArgParseArgument::INTEGER, "INT"));

    addOption(parser, seqan::ArgParseOption(
        "lc", "lc", "Low complexity filter value (percent)",
        seqan::ArgParseArgument::INTEGER, "INT"));

    // Parse command line.
    seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);
    // If parsing was not successful then exit with code 1 if there were errors.
    // Otherwise, exit with code 0 (e.g. help was printed).
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res == seqan::ArgumentParser::PARSE_ERROR;

    // Declare containing variables
    unsigned k = 16;   // default Kmer size
    unsigned lc = 100; // Default low complexity filter (default: none)
    std::string fileName = "example.fa"; // filename
    // Extract option values
    getOptionValue(k, parser, "k");
    getOptionValue(lc, parser, "lc");
    getArgumentValue(fileName, parser, 0);
    // if(lc<25){
    //     std::cout << "Complexity can not be below 25% (nothing would come out...)" << std::endl;
    //     return 1;
    // }
    // counting k-mers
    countKmer(fileName, k, lc);
    
return 0;
}
