
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


typedef std::map< Infix<String<Dna> >::Type , int>  TKmerCounter;




DnaString loadGenome(auto const & filename)
{	
	// Parsing file
	CharString seqFileName = filename;
    StringSet<CharString> ids;
    StringSet<Dna5String> seqs;
    SeqFileIn seqFileIn(toCString(seqFileName));
    readRecords(ids, seqs, seqFileIn);

    // Creating the sequence to index
    DnaString genome;
    for (unsigned i = 0; i < length(ids); ++i){
        append(genome,seqs[i]);
    }
    return(genome);
  }

void countKmer(std::string fileName, int k){

	DnaString genome = loadGenome(fileName);
	// storage
	TKmerCounter kmerSet;
	// temporary storage
	Infix<String<Dna> >::Type inf;

	for( int i = 0; i < length(genome) - k + 1; i++){

		inf = infix(genome, i, i+k);
		if (kmerSet.count(inf) != 0){
			kmerSet[inf] += 1;
		}
		else{
			kmerSet[inf] = 1;
		}
	}

	std::vector<std::pair< Infix<String<Dna> >::Type, int > > kmerVec(std::make_move_iterator(kmerSet.begin()), std::make_move_iterator(kmerSet.end()));

	std::sort(kmerVec.begin(), kmerVec.end(), [](auto x, auto y){ return x.second > y.second;} );

	for (auto elem : kmerVec){
		//std::cout << elem.first << " : " <<elem.second << std::endl;
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
