k#include <iostream>
#include <cstdlib>
#include <math.h>
#include <unordered_map>
#include <chrono>

#include <seqan/index.h>
#include <seqan/seq_io.h>
#include <seqan/sequence.h>
#include <seqan/basic.h>
#include <seqan/stream.h>
#include <seqan/arg_parse.h>


using namespace seqan;


static char bases[] = { 'A', 'C', 'G', 'T' };

// Comparison 
typedef std::unordered_map<std::string, int> resultStore;
typedef std::pair<std::string, int> countResult;
struct CompareSecond
{
    bool operator()(const countResult& left, const countResult& right) const
    {
        return left.second < right.second;
    }
};

// Function to update the result array
countResult updateMin(resultStore &results)
{
    countResult currentMin = *min_element(results.begin(), results.end(), CompareSecond());
    return(currentMin);
}


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

    //Index<DnaString, BidirectionalIndex<FMIndex<void, TFastConfig> > > index(genome);
    typedef FastFMIndexConfig<void, uint32_t, 2, 1> TFastConfig;
    

    // Opening kmerFile
    //CharString kmFileName = kmerFile;

    // Parsing file
    StringSet<CharString> kmIds; // not used, but needed
    //StringSet<DnaString> kmSeqs;
    StringSet<CharString> kmSeqs;
    SeqFileIn kmSeqFileIn(toCString(kmerFile));
    readRecords(kmIds, kmSeqs, kmSeqFileIn);



    // Constants (should be in upper case)
    const int maxErr  =  2;              // Max number of errors, argument of the function
    const int nbRead  = length(ids);     // Number of reads in the file
    const int lenRead = length(seqs[0]);  // size of the sequence to search in
    // Variables used by each thread
    int  limit, th, sum ;
    DnaString seq;
    std::string currentSeq;

    // Shared result storage

    resultStore results;
    std::string minElement;
    int minValue;
    int count[nbRead];

    int range = length(kmSeqs);
    bool isFilled = false;
    std::mutex mut;

    auto delegatePara = [ &mut, &results, &minValue, &minElement, &nbRead, &lenRead, &isFilled, & nbStore](auto & iter, CharString const & needle, uint8_t errors)
    {
        std::lock_guard<std::mutex> lck(mut); // critical section below this line
        // counting presence / absence of k in each sequence.
        unsigned count[nbRead];
        std::fill(count, count+nbRead, 0);

        for (auto occ : getOccurrences(iter))
        {
            count[occ/lenRead] = 1;
        }
        
        // counting total occurence
        int sum = 0;
        for (auto value: count)
        {
            sum += value;
        }

        // actualising results
        
        char* seq = toCString(needle);

        // std::cout << needle << std::endl;
        // for(auto elem: count)
        // {
        //     std::cout << elem << " ";
        // }
        // std::cout << std::endl;

        if(isFilled)
        {
            if(sum > minValue)
            {
                // updating map if required
                results.erase(minElement);
                results[seq] = sum;
                countResult currentMin = updateMin(results);
                minElement = currentMin.first;
                minValue   = currentMin.second;
            }
        }
        else
        {   
            
            results[seq] = sum; // filling map
            // searching and storing min element
            countResult currentMin = updateMin(results);
            minElement = currentMin.first;
            minValue   = currentMin.second;
            if(length(results) == nbStore )
            {
                isFilled = true;
            }
        }
        
    };

    omp_set_num_threads(nbThread);
    #pragma omp paralell shared(range, kmSeqs, mut, results, minValue, minElement, nbRead, lenRead, isFilled) private(delegatePara)
    {   
        // creating one index per thread (avoid parallelisation problem using seqan 4.0)
        Index<DnaString, BidirectionalIndex<FMIndex<TFastConfig> > > index(genome);

       
        std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
        #pragma omp for
        for( int i = 0; i < range; i++ )
        {

            find<0, maxErr>(delegatePara, index, kmSeqs[i] , EditDistance());
            if(i%100000  == 0)
            {
                #pragma omp critical
                {
                    std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
                    auto duration = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();
                    std::cout << i << "/" << range << " : " << duration << " " << omp_get_thread_num() << std::endl;
                }
            }


        }
    }

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