#include <seqan/index.h>
#include <seqan/seq_io.h>
#include <seqan/sequence.h>  
#include <seqan/arg_parse.h>

#include <iostream>
#include <cstdlib>  
#include <fstream>
#include <cmath>
#include <stdlib.h>
#include <time.h>
#include <stdexcept>
#include <unordered_map>



using namespace seqan;

const std::string DNA = "ACGT";

// Setting the index
typedef FastFMIndexConfig<void, uint32_t, 2, 1> TFastConfig;
using index_t = Index<StringSet<DnaString>, BidirectionalIndex<FMIndex<void,TFastConfig> > >;

// counter type, using unordered map.
typedef std::unordered_map<uint,uint> counter;
// pair vector
typedef std::vector<std::pair<unsigned,unsigned> > pair_vector;
// type of sequence set
typedef StringSet<DnaString> sequence_set_type;
// vector of boolean used to keep track of kmer positions count.
typedef std::vector<bool>  bit_field;



// Starting time
const auto boot_time = std::chrono::steady_clock::now();

/**
    Convert a Seqan DnaString to unsigned int.
    SeqAn already store DnaString in 2 bit representation,
    but it is not easy to access it "as is" or to use it
    as key in a hash map.
    @param The sequence to convert
    @return the kmer in 2 bit format, as an unsigned int.
*/
inline unsigned dna2int(DnaString seq){
    
    unsigned value = 0;
    for(auto c : seq){
        value = value << 2 | uint8_t(c);    
    }
    return(value);
}

/**
    Convert an unsigned int back to Seqan DnaString
    @param the integer to convert
    @param the size of the kmer
    @return the kmer in SeqAn DnaString format
*/
inline DnaString int2dna(unsigned value, uint8_t k){
    std::string seq = "";
    for(int i = 0; i< k; i++){
        seq = DNA[value & 3] + seq;
        value >>= 2;

    }
    return(DnaString(seq));
}


/**
    Shortcut to print text in stdout with a time stamp.
    @param whatever you want to print.
*/
template<typename TPrintType>
void print(TPrintType text)
{
    const auto milis = std::chrono::duration <double, std::milli>(std::chrono::steady_clock::now() - boot_time).count();
    std::cout << "[" << milis << " ms]\t" << text << std::endl;
}


/**
    Shortcut to print list of kmer counts (using either pairs or map )
    @param the counter to export (pair or map of kmer associated to their count)
    @param k, the size of the kmer counted (needed for conversion 2bit representation back to DNA)
*/
template<typename TPrintType>
void printCounters(TPrintType & pvec, uint8_t k){
    for(auto it = pvec.begin(); it!= pvec.end(); ++it)
    { 
        std::cout << int2dna(it->first,k) << " " << it->second << "\n";
    }   
}

/**
    Export a counter to a file
    @param the counter to export (pair or map of kmer associated to their count)
    @param k, the size of the kmer counted (needed for conversion 2bit representation back to DNA)
    @param the path to the outputfile
*/
template<typename TPrintType>
void exportCounter(TPrintType & pvec, uint8_t k, std::string output){

    std::ofstream outputFile;
    outputFile.open (output);
    if(outputFile.is_open()){
        for(auto it = pvec.begin(); it!= pvec.end(); ++it)
        { 
            outputFile << int2dna(it->first,k) << "\t" << it->second << "\n";
        }   
        outputFile.close();
    }
    else{
        print("COULD NOT OPEN FILE " + output);
    }
}

/**
    Adjust low complexity threshold value to kmer size
    @param low complexity threshold for a size of kmer
    @param kmer size for this threshold
    @param new k-mer size
    @return the appropriate threshold value for the new kmer size.
*/
float adjust_threshold(float c_old, uint8_t k_old, uint8_t k_new ){
    float c_new = c_old * float( std::pow(k_new - 2 + 1,2) /  std::pow(k_old - 2 + 1,2));
    return(c_new);
}

/**
    Compute the sum of a vector (of int for example).
    @param the vector
    @return the sum
*/
template<typename TVector>
inline unsigned vectorSum(TVector vec)
{
    unsigned res = 0;
    for(auto it: vec)
    {
        res += it;
    }
    return(res);
}

/**
    Check the complexity of a kmer
    @param the kmer to test, in 2 bit representation, cast as an unsigned int
    @param k the size of said kmer
    @param threshold for the low complexity filter.
    @return True if the kmer contains low complexity region
*/
inline  bool haveLowComplexity(unsigned kmer, uint8_t k, float threshold){
    
    unsigned counts[16] = { 0 }; // 16 possibles dimers
    // reading using sliding window of 2
    for(int i = 0; i < k-1; i++){
        // storing dimers as 2 * 2 bits
        uint8_t c =  kmer & 16;
        // removing last element of the k-mer
        kmer >>=2;
        // updating value of dimer count
        counts[c]++;
    }
    float s = 0;
    size_t sum = 0;
    for(auto v:counts){
        sum +=  v * (v-1);  
    }
    s =  sum / float(2 * (k-2));
    return s>= threshold;
}

/**
    Return the first top kmers, ranked by count, descending.
    @param kmer count map (see count_kmer and error count)
    @param number of kmers to return
    @return vector of pair containing the most frequent kmers and their associated count.
*/
pair_vector get_most_frequent(counter & count_map, unsigned limit){

    pair_vector kmer_vec(std::make_move_iterator(count_map.begin()), std::make_move_iterator(count_map.end()));
    std::sort(kmer_vec.begin(), kmer_vec.end(), [](auto x, auto y){ return x.second > y.second;} );
    if(kmer_vec.size() > limit){
        kmer_vec.resize(limit);
    }
    return(kmer_vec);
}


/**
    Sample the sequences set and return requested samples cut to size
    @param Set of sequences to sample (SeqAn StringSet of DnaString)
    @param Number of sequences to sample
    @param Size of the sampled sequence.
    @return A set of sample sequences cut to size.
*/
sequence_set_type sampleSequences(sequence_set_type & sequence_set, unsigned sample_size, unsigned cut_size){
    sequence_set_type sample;
    
    // Initialising the random seed
    srand(time(0));
    
    // Checking if we can sample the requested number of sequences, else return the whole set
    unsigned nb_sample = sample_size;
    unsigned sequence_set_size = length(sequence_set);
    if(sample_size > sequence_set_size){
        std::cout << "Sequence set too small for the requested sample size\n";
        std::cout << "The whole set will be used." ;
        nb_sample = sequence_set_size;
    }

    // Building vector with all possible seq indice
    std::vector<int> v(sequence_set_size) ; 
    std::iota(std::begin(v), std::end(v), 0);

    // and applying random shuffling to said vector.
    std::random_device rd;
    std::mt19937 g(rd());
    std::shuffle(v.begin(), v.end(), g);

    // counters
    unsigned nb_seq = 0;
    unsigned i = 0;
    unsigned seq_id;

    // Fetching the random sequences
    while(nb_seq < nb_sample and i < sequence_set_size ){
        
        // current sequence id
        seq_id = v[i];

        // we also need to check that the sequence is
        // at least as long as the requested sample length
        if( length(sequence_set[ seq_id ]) >= cut_size ){
            appendValue(sample, prefix(sequence_set[ seq_id ],cut_size));
            nb_seq +=1;
        }
        i++;
    }
    print("Sampled " + std::to_string(length(sample)) + " sequences");
    return(sample);
}


/**
    Perform a simple exact count of all the k-mers from a sample set of sequences
    @param Set of sequences (SeqAn StringSet of DnaString)
    @param k, size of the kmers
    @param threshold of the low complexity filter
    @return a map of the kmer count, with a kmer hash as key.

*/
counter count_kmers(sequence_set_type & sequences, uint8_t k, float threshold){

    counter count;
    unsigned base = std::pow(2,(2*k))-1;
    for(auto seq: sequences){
        // First kmer
        unsigned n = dna2int(DnaString(prefix(seq,k-1)));
        int i = k-1;
        while(i < length(seq)){ 
            
            n <<= 2; 
            n = (n & base) |  size_t(seq[i]);
            if(not haveLowComplexity(n,k,threshold)){
                count[n] +=1 ;
            }           
            i++;
        }
    }
    return(count);
}


/**
    Search and count a list of kmer in a set of sequences, at a Levenstein distance of at most 2.
    @param Set of sequences (SeqAn StringSet of DnaString)
    @param the previous count of exact kmer (the kmer list)
    @param number of thread to use
    @param k, size of the kmers
    @return a map of the kmer count, with a kmer hash as key.

*/
counter errorCount( sequence_set_type & sequences, pair_vector & exact_count, uint8_t nb_thread, uint8_t k){
    
    const uint8_t MAXERR = 2; // Max number of errors, need to be fixed at compile time

    unsigned sample_size = length(sequences);

    print("Preparing index");
    index_t  index(sequences);
    print("Creating index");
    indexCreate(index);
    
    // Result storage
    counter results;

    // setting number of parallel thread
    omp_set_num_threads(nb_thread);
    print("Starting approximate counting");
    #pragma omp parallel shared(index, results)
    {
        // local variable to keep track of kmer occurences
        std::array<bit_field,3> tcount;

        // Delegate function for SeqAn find function (process occurences)
        auto delegateParallel = [& tcount](auto & iter, const DnaString & needle, int errors)
        {
            for (auto occ : getOccurrences(iter)){
                
                unsigned read_id = getValueI1(occ);
                tcount[errors][read_id] = true;
                }
        };

        #pragma omp for schedule(dynamic)
        for(int km_id=0; km_id<length(exact_count); km_id++)
        {
            
            /*  Supposing a kmer is very unlikely to be twice in the same read start
                we just store the read id in which kmer have been found using a bit field
                which size is the number of reads.
                Since we search at 2 error, we can either:
                    - Find the kmer exactly
                    - Find it at 1 error
                    - Find it at 2 errors
                So there is one bit field per number of error.
            */
            for(int i=0; i<3; i++){
                tcount[i] = bit_field(sample_size,false);
            }
            // 2 bit encoded kmer as unsigned int
            unsigned kmer = exact_count[km_id].first;
            // ressearch, filling tcount
            find<0, MAXERR >(delegateParallel, index, int2dna(kmer,k), EditDistance() );

            
            // computing total number of occurences
            unsigned total = 0;
            for(auto bit_count: tcount){
                total +=  vectorSum(bit_count);
            }
            // Updating global counter
            #pragma omp critical
            results[kmer] = total;
                
        }
    }
    return(results);
}




int main(int argc, char const ** argv)
{

    // Setup ArgumentParser.
    seqan::ArgumentParser parser("adaptFinder");

    addArgument(parser, seqan::ArgParseArgument(
        seqan::ArgParseArgument::STRING, "input filename"));

    addOption(parser, seqan::ArgParseOption(
        "lc", "low_complexity", "low complexity filter threshold (for k=16), default 1.5",
        seqan::ArgParseArgument::DOUBLE, "kmer filename"));
    
    addOption(parser, seqan::ArgParseOption(
        "sn", "sample_n", "sample n sequences from dataset, default 10k sequences",
        seqan::ArgParseArgument::INTEGER, "INT"));

    addOption(parser, seqan::ArgParseOption(
        "sl", "sample_length", "size of the sampled portion, default 100 bases",
        seqan::ArgParseArgument::INTEGER, "INT"));

    addOption(parser, seqan::ArgParseOption(
        "nt", "nb_thread", "Number of thread to work with, default is 4",
        seqan::ArgParseArgument::INTEGER, "INT"));

    addOption(parser, seqan::ArgParseOption(
        "k", "kmer_size", "Size of the kmers, default is 16",
        seqan::ArgParseArgument::INTEGER, "INT"));

    addOption(parser, seqan::ArgParseOption(
        "lim", "limit", "limit the number of kmer used after initial counting, default is 500",
        seqan::ArgParseArgument::INTEGER, "INT"));

     addOption(parser, seqan::ArgParseOption(
        "v", "verbosity", "Level of details printed out",
        seqan::ArgParseArgument::INTEGER, "INT"));

    addOption(parser, seqan::ArgParseOption(
        "e", "exact_file", "path to export the exact k-mer count, if needed. Default: no export",
        seqan::ArgParseArgument::STRING, "exact count output file"));

    addOption(parser, seqan::ArgParseOption(
        "o", "out_file", "path to the output file, default is ./out.txt",
        seqan::ArgParseArgument::STRING, "output file"));

    // Parser command line.
    seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);

    // If parsing was not successful then exit with code 1. if there were errors.
    // Otherwise, exit with code 0 (e.g. help was printed).
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res == seqan::ArgumentParser::PARSE_ERROR;

    // Extract option values and print them.
    std::string output = "out.txt";     // output file
    std::string exact_out;   // exact count output file
    unsigned nb_thread = 4;  // default number of thread
    unsigned k = 16;         // kmer size, 2<= k <= 32
    unsigned sl = 100 ;      // sequence sampling size
    unsigned sn = 10000;     // number of sequence sampled
    unsigned limit = 500;    // number of kmers to keep.
    double lc = 1.5;         // low complexity filter threshold, allow all known adapters to pass.
    unsigned v = 0;

    // Fetching values
    getOptionValue(limit, parser, "lim");
    getOptionValue(lc, parser, "lc");
    getOptionValue(k, parser, "k");
    getOptionValue(v, parser, "v");
    getOptionValue(sl, parser, "sl");
    getOptionValue(sn, parser, "sn");
    getOptionValue(nb_thread, parser, "nt");
    getOptionValue(output, parser, "o");
    getOptionValue(exact_out, parser, "e");
    std::string input_file;
    getArgumentValue(input_file, parser, 0);


    // checking value for k
    if( k<2 or k>32 ){
        throw std::invalid_argument("kmer size must be between 2 and 32 (included)");
    }
    
    // adjusting low complexity to kmer size
    lc = adjust_threshold( lc, 16, k );
    print("LC filter adjusted to " + std::to_string(lc));
    // Parsing input fasta file
    StringSet<CharString> ids;
    StringSet<DnaString> seqs;
    print("Parsing FASTA file");
    SeqFileIn seqFileIn(toCString(input_file));
    readRecords(ids, seqs, seqFileIn);

    // sample and cut sequences to required length
    print("Sampling");
    sequence_set_type sample = sampleSequences(seqs, sn, sl);
    
    // counting k-mers on the sampled sequences
    print("Exact k-mer count");
    counter count = count_kmers(sample, k, lc);
    
    // keeping most frequents kmers
    print("Number of kmer found / kept");
    print(count.size());
    pair_vector first_n_vector = get_most_frequent(count, limit);
    print(count.size());

    // Exporting exact kmer count, if required
    if(not exact_out.empty() ){
        print("Exporting exact kmer count");
        exportCounter(first_n_vector, k, exact_out);
        //exportCounter(get_most_frequent(count, limit),k,exact_out);   
    }

    // Counting with at most 2 errors
    print("Approximate k-mer count");
    counter error_counter = errorCount(sample, first_n_vector, nb_thread, k);
    pair_vector sorted_error_count = get_most_frequent(error_counter,limit);

    //printCounters(sorted_error_count,k);
    print("Exporting approximate count");
    exportCounter(sorted_error_count,k,output);
    print("Done");
    return 0;
}