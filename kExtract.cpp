
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


typedef std::map< Dna5String, int>  TKmerCounter;




StringSet<DnaString> loadGenome(auto const & filename)
{	
	// Parsing file
	CharString seqFileName = filename;
    StringSet<CharString> ids;
    StringSet<Dna5String,Owner<> > seqs;
    SeqFileIn seqFileIn(toCString(seqFileName));
    readRecords(ids, seqs, seqFileIn);

    return(seqs);
  }

void countKmer(std::string fileName, int k){

	StringSet<Dna5String> sequences = loadGenome(fileName);
	// storage
	TKmerCounter kmerSet;

    for(int j = 0; j< length(sequences); j++)
    {
        Dna5String seq = sequences[j];
    	for( int i = 0; i < length(seq) - k + 1; i++){

    		Infix<Dna5String>::Type inf = infix(seq, i, i+k);
            Dna5String km = Dna5String(inf);

            if ( kmerSet.count(km) !=0 ){
    			kmerSet[km] += 1;
            }
    		
    		else{
    			kmerSet[km] = 1;
                //testVec.push_back(inf);
    		}
    	}
    }
	std::vector<std::pair< Dna5String, int > > kmerVec(std::make_move_iterator(kmerSet.begin()), std::make_move_iterator(kmerSet.end()));
    std::sort(kmerVec.begin(), kmerVec.end(), [](auto x, auto y){ return x.second > y.second;} );


	for (auto elem : kmerVec)
    {
		std::cout << elem.first << "\t" <<elem.second << std::endl;
		//std::cout << ">_" << elem.first << "_" << elem.second  << std::endl;
        //std::cout << elem.first << std::endl;
	// }
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

    // Parse command line.
    seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);
    // If parsing was not successful then exit with code 1 if there were errors.
    // Otherwise, exit with code 0 (e.g. help was printed).
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res == seqan::ArgumentParser::PARSE_ERROR;

    // Extract option values and print them.
    unsigned k = 16;   // default Kmer size
    getOptionValue(k, parser, "k");

    std::string fileName = "example.fa";
    getArgumentValue(fileName, parser, 0);
    countKmer(fileName, k);

return 0;
}

// pseudocode

// for fragment in genome:
// 	if fragment in ( list )
// 		list[fragment].count += 1
// 	else:
// 		list[fragment].count = 1

// To sort stuff by value:
// std::sort(items.begin(), items.end(), cmp);
