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



using namespace seqan;

const std::string DNA = "ACGT";

// counter type, using unordered map.
typedef std::unordered_map<uint,uint> counter;
// Setting the index
typedef FastFMIndexConfig<void, uint32_t, 2, 1> TFastConfig;
// pair vector
typedef std::vector<std::pair<unsigned,unsigned> > pair_vector;
// type of sequence set
typedef StringSet<DnaString> sequenceSetType;
// vector of boolean used to keep track of kmer positions count.
typedef std::vector<bool>  bitField;

using index_t = Index<StringSet<DnaString>, BidirectionalIndex<FMIndex<void,TFastConfig> > >;


inline unsigned dna2int(DnaString seq){
    
    unsigned value = 0;
    for(auto c : seq){
        value = value << 2 | uint8_t(c);    
    }
    return(value);
}

inline DnaString int2dna(unsigned value, uint8_t k){
    std::string seq = "";
    for(int i = 0; i< k; i++){
        seq = DNA[value & 3] + seq;
        value >>= 2;

    }
    return(DnaString(seq));
}


const auto boot_time = std::chrono::steady_clock::now();
// Shortcut to print text in stdout
template<typename TPrintType>
void print(TPrintType text)
{
    const auto milis = std::chrono::duration <double, std::milli>(std::chrono::steady_clock::now() - boot_time).count();
    std::cout << "[" << milis << " ms]\t" << text << std::endl;
}


// Shortcut to print list of kmer counts (using either pairs or map )
template<typename TPrintType>
void printCounters(TPrintType & pvec, uint8_t k){
    for(auto it = pvec.begin(); it!= pvec.end(); ++it)
    { 
        std::cout << int2dna(it->first,k) << " " << it->second << "\n";
    }   
}

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

pair_vector get_most_frequent(counter & count_map, unsigned limit){

    pair_vector kmer_vec(std::make_move_iterator(count_map.begin()), std::make_move_iterator(count_map.end()));
    std::sort(kmer_vec.begin(), kmer_vec.end(), [](auto x, auto y){ return x.second > y.second;} );
    if(kmer_vec.size() > limit){
        kmer_vec.resize(limit);
    }
    return(kmer_vec);
}


counter count_kmers(sequenceSetType & sequences, uint8_t k, float threshold){

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


// Cut all sequences from the set to a precise length
sequenceSetType cutSequences(sequenceSetType & sequence_set, unsigned cut_size ){
    sequenceSetType cutSequences;
    for(auto seq: sequence_set){
        unsigned l = length(seq);
        appendValue(cutSequences, prefix(seq,std::min(cut_size,l) ) );
    }
    return(cutSequences);
}

// Select sample_size sequences and return them
sequenceSetType sampleSequences(sequenceSetType & sequence_set, unsigned sample_size){
    sequenceSetType sample;
    
    // INIT THE RANDOM SEED
    srand(time(0));
    
    // simple checksum to check if the random seeding works properly.
    unsigned checksum = 0; 

    // Input parameters and control variables.
    unsigned seq_length = length(sequence_set);
    unsigned nb_sample = std::min(sample_size,seq_length);
    unsigned seq_id;

    // Building vector with all possible seq indice
    std::vector<int> v(seq_length) ; 
    std::iota (std::begin(v), std::end(v), 0);

    // and applying random shuffling to said vector.
    std::random_device rd;
    std::mt19937 g(rd());
 
    std::shuffle(v.begin(), v.end(), g);

    // Finally reporting the sequences of interest.
    for(int i = 0; i< nb_sample; i++){
        seq_id = v[i];
        DnaString seq = sequence_set[ seq_id ];
        appendValue(sample,seq);
        checksum += seq_id;
    }
    std::cout << "Checksum: " << checksum << "\n";
    std::cout << "Exporting " << length(sample) << " sequences.\n";
    return(sample);
}


// Search and count a kmer list in a fasta file, at a levenstein distance of at most  2.
counter errorCount( sequenceSetType & sequences, pair_vector & exact_count, unsigned sample_size, uint8_t nb_thread, uint8_t k){
    
        const uint8_t MAXERR = 2; // Max number of errors

        print("PREPARING INDEX");
        index_t  index(sequences);
        print("CREATING INDEX");
        indexCreate(index);
        
        // Result storage
        counter results;

        // setting number of parallel thread
        omp_set_num_threads(nb_thread);
        print("STARTING APPROXIMATE COUNTING");
        #pragma omp parallel shared(index)
        {
            // local variable to keep of kmer occurences
            std::array<bitField,3> tcount;
            auto delegateParallel = [& tcount](auto & iter, const DnaString & needle, int errors)
            {
                for (auto occ : getOccurrences(iter)){
                    
                    unsigned read_id = getValueI1(occ);
                    tcount[errors][read_id] = true;
                    }
            };

        #pragma omp for schedule(dynamic)
        for(int km=0; km<length(exact_count); km++)
        {
            
            
            // Supposing a kmer is very unlikely to be twice in the same read start
            // we just store the read id in which kmer have been found using a bit field.
            // the same size as the number of reads
            for(int i=0; i<3; i++){
                tcount[i] = bitField(sample_size,false);
            }
            // 2 bit encoded kmer as unsigned int
            unsigned kmer = exact_count[km].first;
            // ressearch, filling tcount
            find<0, MAXERR >(delegateParallel, index, int2dna(kmer,k), EditDistance() );

            
            // computing total number of OC
            unsigned total = 0;
            for(auto bit_count: tcount){
                total +=  vectorSum(bit_count);
            }
            // critical section, updating counter
            #pragma omp critical
            results[kmer] = total;
                
        }
    }
    print("DONE");
    return(results);
}




int main(int argc, char const ** argv)
{


    // Setup ArgumentParser.
    seqan::ArgumentParser parser("polo_counter");

    addArgument(parser, seqan::ArgParseArgument(
        seqan::ArgParseArgument::STRING, "input filename"));

    addOption(parser, seqan::ArgParseOption(
        "lc", "low_complexity", "low complexity filter threshold (for k16), default 1.5",
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

    // Parse command line.
    seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);

    // If parsing was not successful then exit with code 1 if there were errors.
    // Otherwise, exit with code 0 (e.g. help was printed).
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res == seqan::ArgumentParser::PARSE_ERROR;

    // Extract option values and print them.
    std::string output = "out.txt";     // output file
    std::string exact_out;     // exact count output file
    unsigned nb_thread = 4;  // default number of thread
    unsigned k = 16;        // kmer size
    unsigned sl = 100 ;     // sequence sampling size
    unsigned sn = 10000;    // number of sequence sampled
    unsigned limit = 500;   // number of kmers to keep.
    double lc = 1.5;        // low complexity filter threshold, allow all known adapters to pass.
    unsigned v = 0;
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
    
    
    // Parsing input fasta file
    StringSet<CharString> ids;
    StringSet<DnaString> seqs;
    print("PARSING FILES");
    SeqFileIn seqFileIn(toCString(input_file));
    readRecords(ids, seqs, seqFileIn);


    print("SAMPLING");
    // cut sequences to required length
    sequenceSetType sample = cutSequences(sampleSequences(seqs,sn),sl);
    print("EXACT K-MER COUNT");
    // counting k-mers on the sampled sequences
    counter count = count_kmers(sample, k, lc);
    // keeping most frequents
    print(count.size());
    pair_vector first_n_vector = get_most_frequent(count, count.size());
    
    if(not exact_out.empty() ){
        print("EXPORTING EXACT COUNT");
        //exportCounter(first_n_vector,k,exact_out);
        exportCounter(get_most_frequent(count, limit),k,exact_out);
        
    }

    print("APPROXIMATE K-MER COUNT");
    // Counting at 2 errors
    counter error_counter = errorCount(sample,first_n_vector, sn, nb_thread, k);
    pair_vector sorted_error_count = get_most_frequent(error_counter,limit);
    //printCounters(sorted_error_count,k);
    print("EXPORTING APPROXIMATE COUNT");
    exportCounter(sorted_error_count,k,output);
    
    return 0;
}