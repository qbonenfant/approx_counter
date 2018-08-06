#include <iostream>
#include <cstdlib>

#include <seqan/index.h>
#include <seqan/seq_io.h>
#include <seqan/sequence.h>
#include <seqan/basic.h>
#include <seqan/stream.h>
#include <seqan/arg_parse.h>


using namespace seqan;

// Defining counter type
// I use an unordered map, using Dna string as keys, which is associated to an array of three vectors
// Each vector  are 'number of read long', and contain a bool defining if yes or no, the key is found
// in the n th sequence, at a given error rate( here, either 0, 1 or 2 )
typedef std::unordered_map<DnaString, std::array<std::vector<bool>, 3> > counter;


// subtype for elements
typedef std::pair<DnaString, std::array<std::vector<bool>, 3> > counterElem;

// Setting the index
typedef FastFMIndexConfig<void, uint32_t, 2, 1> TFastConfig;




// Sale, mais obligé vu que SeqAn ne le fait pas tout seul //
// Defining a hash function for DnaString
namespace std {
    template<>
    class hash<DnaString> {
    public:
        size_t operator() (const DnaString& v) const 
        {
            // TODO: Voir si on peut optimiser...
            size_t result = 0;
            for(auto c : v) {
                result = size_t(c) + 31*result;
            }
            return result;
        }
    };
}

// Shortcut to print text in stdout
template<typename TPrintType>
void print(TPrintType text)
{
    std::cout << text << std::endl;
}
&

template<typename TVector>
int vectorSum(TVector vec)
{
    int res = 0;
    for(auto it: vec)
    {
        res += it;
    }
    return(res);
}


// Function used to print results
void printCounter(counter count)
{
    std::vector<std::pair< std::array<int,3>, DnaString > > toPrint;
    // converting
    for(auto it = count.begin(); it!= count.end(); ++it)
    {   
        std::array<int,3> tmp = {0,0,0};
        for(int i = 0; i<3 ; i++)
        {
            tmp[i] = vectorSum(it->second[i]);
        }

        toPrint.push_back( std::make_pair(tmp, it->first ) );
    
    }
    // sorting by count
    std::sort(toPrint.begin(),toPrint.end(), [](auto e1, auto e2){ 
        int v1 = vectorSum(e1.first);
        int v2 = vectorSum(e2.first);
        return v1>v2;} );
    

    for(auto it = toPrint.begin(); it!= toPrint.end(); ++it)
    {
        for(auto e: it->first)
        {
            std::cout << e << "\t";
        }
        std::cout << it->second << std::endl;
    }
}




// template< typename TElement1, typename TElement2 >
// bool minCounter(TElement1 e1, TElement2 e2 )
// {
//     for(int i = 0; i< 3; i++)
//     {
//         int v1 = vectorSum(e1.second[i]);
//         int v2 = vectorSum(e2.second[i]);
//         if(v1!=v2)
//         {
//             return(v1<v2);
//         }
//     }
// }


// Search and count a kmer list in a fasta file, at at most a levenstein distance of 2.
void findAdapt(const std::string& filename, std::string& kmerFile, const int& nbErr, const int& nbStore, const int& nbThread ){
    
     // Opening input fasta file

    // Parsing file
    StringSet<CharString> ids;
    StringSet<DnaString> seqs;
    SeqFileIn seqFileIn(toCString(filename));
    readRecords(ids, seqs, seqFileIn);


    // Creating the sequence to index
    DnaString genome;
    for (unsigned i = 0; i < length(ids); ++i){
        append(genome,seqs[i]);
    }

    

    // Opening kmerFile

    // Parsing file
    StringSet<CharString> kmIds; // not used, but needed
    StringSet<CharString> kmSeqs;
    SeqFileIn kmSeqFileIn(toCString(kmerFile));
    readRecords(kmIds, kmSeqs, kmSeqFileIn);


    // Constants (should be in upper case)
    const int maxErr  = 2;               // Max number of errors, argument of the function
    const int nbRead  = length(ids);     // Number of reads in the file
    const int lenRead = length(seqs[0]); // size of the sequence to search in

    counter results;
    bool count[nbRead];
    // dummy sequence, used for research initialisation
    DnaString dummy = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";



    // mutex for atomic operations
    std::mutex mtx;
    auto delegateParallel = [&mtx, &count, &results, &nbRead, &lenRead, &nbStore](auto & iter, const DnaString & needle, int errors)
    {

        std::lock_guard<std::mutex> lck(mtx); // critical section below this line
        std::fill(count, count+nbRead, false);
        
        for (auto occ : getOccurrences(iter)){
            // Testing if the occurence is on one read only
            // (and not between two reads)
            if(occ%lenRead + length(needle) - errors <= lenRead )
            {
                count[occ/lenRead] = true;
                //std::cout << omp_get_thread_num() << " " << occ << " " << needle << " " << errors << std::endl;
            }
            // // For debug only
            // else
            // {
            //     print("Read Overlap");
            // }
        }
        

        // If new element
        if( results.find(needle) == results.end()  )
        {
            // Check if we already have enough elemnts
            if( results.size() >= nbStore +1 )
            {
            
                // auto minIt = min_element(results.begin(), results.end(),
                // [](decltype(results)::value_type& l, decltype(results)::value_type& r) -> bool { return l.second < r.second; });
                auto minIt = min_element(results.begin(), results.end(),
                [](decltype(results)::value_type& l, decltype(results)::value_type& r) -> bool {
                    int v1 = 0;
                    int v2 = 0;
                    for(int i = 0; i< 3; i++)
                    {
                        v1 += vectorSum(l.second[i]);
                        v2 += vectorSum(r.second[i]);
                    return(v1<v2);
                    }
                });

                results.erase(minIt);

            }
            // add new element so it can be filled
            for(int i =0; i<3 ; i++)
            {   
                results[needle][i] =  std::vector<bool>(nbRead,false);
            }
            
            
        }
         
        // filling current element
        for(int i = 0; i < nbRead ; i++ )
        {

            results[needle][errors][i] = count[i] | results[needle][errors][i] ;
        }
    };


    Index<DnaString, BidirectionalIndex<FMIndex<TFastConfig> > > index(genome);
    // trying to perform a fake ressearch to "fix" the problem with seqan find parallel
    // using long homopolymer, that won't modify results.
    
    find<0, maxErr>(delegateParallel, index, dummy , EditDistance() );
    //std::cout<<"Dummy search done" << std::endl;

    results.clear();   // clearing potential parasite results

    //std::cout << "\n  DUMMY SEARCH DONE \n" << std::endl;


    omp_set_num_threads(nbThread);    
    find<0, maxErr>(delegateParallel, index, kmSeqs , EditDistance(),Parallel());

    if(results.size() > nbStore)
    {
        auto minIt = min_element(results.begin(), results.end(),
        [](decltype(results)::value_type& l, decltype(results)::value_type& r) -> bool {
            int v1 = 0;
            int v2 = 0;
            for(int i = 0; i< 3; i++)
            {
                v1 += vectorSum(l.second[i]);
                v2 += vectorSum(r.second[i]);
            return(v1<v2);
            }
        });

        results.erase(minIt);
    }
    printCounter(results);
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

    addOption(parser, seqan::ArgParseOption(
        "ne", "nbErr", "Number of authorized errors (/!\\not used)",
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
    unsigned nbErr = 2;     // default number of errors
    getOptionValue(k, parser, "kf");
    getOptionValue(nbThread, parser, "nt");
    getOptionValue(nbStore,  parser, "ns");
    getOptionValue(nbErr,  parser, "ne");

    std::string text;
    getArgumentValue(text, parser, 0);
    std::cout << text << " " << k << " " <<  nbStore << " " << nbThread << std::endl;
    findAdapt(text,k,nbErr,nbStore,nbThread);

    return 0;
}