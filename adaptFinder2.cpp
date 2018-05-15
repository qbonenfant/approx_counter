#include <iostream>
#include <cstdlib>
#include <chrono>

#include <seqan/index.h>
#include <seqan/seq_io.h>
#include <seqan/sequence.h>
#include <seqan/basic.h>
#include <seqan/stream.h>
#include <seqan/arg_parse.h>


using namespace seqan;


void findAdapt(const std::string& filename, std::string& kmerFile , const int& nbStore, const int& nbThread ){
    
     // Opening input fasta file

    // CharString seqFileName = "/home/cube/Documents/Code/training/SeqAn/adaptateurs/start20k_Filtered50.fa";
    // CharString seqFileName = "/home/cube/Documents/Code/training/SeqAn/adaptateurs/example.fa";
    CharString seqFileName = filename;

    // Parsing file
    StringSet<CharString> ids;
    StringSet<DnaString> seqs;
    SeqFileIn seqFileIn(toCString(seqFileName));
    readRecords(ids, seqs, seqFileIn);


    // Creating the sequence to index
    DnaString genome;
    for (unsigned i = 0; i < length(ids); ++i){
        append(genome,seqs[i]);
    }

    // Setting FMindex


    

    // Opening kmerFile
    //CharString kmFileName = kmerFile;

    // Parsing file
    StringSet<CharString> kmIds; // not used, but needed
    //StringSet<DnaString> kmSeqs;
    StringSet<CharString> kmSeqs;
    SeqFileIn kmSeqFileIn(toCString(kmerFile));
    readRecords(kmIds, kmSeqs, kmSeqFileIn);



    // Constants (should be in upper case)
    const int maxErr  =  1;              // Max number of errors, argument of the function
    const int nbRead  = length(ids);     // Number of reads in the file
    const int lenRead = length(seqs[0]);  // size of the sequence to search in

    int count[nbRead];
    // Defining types to store results
    typedef std::pair<int,DnaString> resultCount;
    typedef std::vector<resultCount> resultStore;
    // Storing counts and sequences
    resultStore results;



    // Setting the index
    //Index<DnaString, BidirectionalIndex<FMIndex<void, TFastConfig> > > index(genome);
    typedef FastFMIndexConfig<void, uint32_t, 2, 1> TFastConfig;


    std::mutex mtx;
    auto delegateParallel = [&count, &results, &mtx, &nbRead,  &lenRead, &nbStore](auto & iter, const DnaString & needle, uint8_t errors)
    {

        std::lock_guard<std::mutex> lck(mtx); // critical section below this line
        std::fill(count, count+nbRead, 0);
        
        for (auto occ : getOccurrences(iter)){
            count[occ/lenRead] = 1;
        }

        int sum = 0;
        for (int i = 0; i<nbRead; i++){
                sum += count[i];
            }
        resultCount currentCount = std::make_pair( sum ,needle);
        
        if(std::find(results.begin(), results.end(), currentCount ) == results.end())
        {
            results.push_back(currentCount);
        }
        if(length(results) > nbStore )
        {
            auto minElem = min_element(results.begin(), results.end());
            results.erase(minElem);
        }

    };


    Index<DnaString, BidirectionalIndex<FMIndex<TFastConfig> > > index(genome);
    // trying to perform a fake ressearch to "fix" the problem with sean find parallel
    // using long homopolymer, that won't modify result, hopefully
    find<0, maxErr>(delegateParallel, index, kmSeqs[0] , EditDistance() );


    omp_set_num_threads(nbThread);    
    find<0, maxErr>(delegateParallel, index, kmSeqs , EditDistance(),Parallel());

    //     if(i%100000  == 0)
    //     {
    //         #pragma omp critical
    //         {
    //             std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
    //             auto duration = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();
    //             std::cout << i << "/" << range << " : " << duration << " " << omp_get_thread_num() << std::endl;
    //         }
    //     }
    

    std::sort(results.begin(),results.end(), [](auto e1, auto e2){ return e1.first> e2.first;} );
    for (auto elem: results)
    {
        std::cout << elem.first << " ; " << elem.second << std::endl;
    }


}

int main(int argc, char const ** argv)
{
    // Setup ArgumentParser.
    seqan::ArgumentParser parser("findAdapt");

    addArgument(parser, seqan::ArgParseArgument(
        seqan::ArgParseArgument::STRING, "input filename"));

    addOption(parser, seqan::ArgParseOption(
        "kf", "kmerFile", "path to the kmer file",
        seqan::ArgParseArgument::STRING, "kmer filename"));
    addOption(parser, seqan::ArgParseOption(
        "nt", "nbThread", "Number of thread to work with",
        seqan::ArgParseArgument::INTEGER, "INT"));

    addOption(parser, seqan::ArgParseOption(
        "ns", "nbStore", "Number of kmer to keep in result list",
        seqan::ArgParseArgument::INTEGER, "INT"));

    // Parse command line.
    seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);

    // If parsing was not successful then exit with code 1 if there were errors.
    // Otherwise, exit with code 0 (e.g. help was printed).
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res == seqan::ArgumentParser::PARSE_ERROR;

    // Extract option values and print them.
    std::string k = "kexample.fa"; // default Kmer file
    unsigned nbStore = 10;  // default report size
    unsigned nbThread = 4;  // default number of thread
    getOptionValue(k, parser, "kf");
    getOptionValue(nbThread, parser, "nt");
    getOptionValue(nbStore,  parser, "ns");

    std::string text;
    getArgumentValue(text, parser, 0);
    std::cout << text << " " << k << " " <<  nbStore << " " << nbThread << std::endl;
    findAdapt(text,k,nbStore,nbThread);

    return 0;
}