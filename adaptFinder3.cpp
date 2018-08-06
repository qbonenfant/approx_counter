
#include <seqan/index.h>
#include <seqan/seq_io.h>
#include <seqan/sequence.h>
#include <seqan/basic.h>
#include <seqan/stream.h>
#include <seqan/arg_parse.h>

#include <iostream>
#include <cstdlib>
#include <fstream>


using namespace seqan;


// Defining counter type
// I use an unordered map, using c strings as keys, which is associated to an array of three vectors
// Each vector  are 'number of read' long, and contain a bool defining if yes or no, the key is found
// in the n th sequence, at a given error rate( here, either 0, 1 or 2 )
typedef std::unordered_map<std::string, std::array<std::vector<bool>, 3> > counter;



// Setting the index
typedef FastFMIndexConfig<void, uint32_t, 2, 1> TFastConfig;



// Shortcut to print text in stdout
template<typename TPrintType>
void print(TPrintType text)
{
    std::cout << text << std::endl;
}


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


void printVec(std::vector<bool> vec){
    for(auto it: vec){
        std::cout << it << "\t";
    }
    std::cout << std::endl;
}

// Function used to sort and print results in file
void printCounter(counter count, std::ofstream &outputFile)
{
    std::vector<std::pair< std::array<int,3>, std::string > > toPrint;
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
            outputFile << e << "\t";
        }
        outputFile << it->second << std::endl;
    }
}

// // Reset the temporary counters
// void tempCountReset(tempCount &tc, int size)
// {   
//     for(int i =0; i<3; i++)
//     {
//         std::fill(tc[i], tc[i]+size, false);
//     }
// }


// Search and count a kmer list in a fasta file, at at most a levenstein distance of 2.
void findAdapt(const std::string& filename, std::string& kmerFile, const int& nbThread, std::ofstream &outputFile ){
    
    print("PARSING FILES");
    // Parsing input fasta file
    StringSet<CharString> ids;
    StringSet<DnaString> seqs;
    SeqFileIn seqFileIn(toCString(filename));
    readRecords(ids, seqs, seqFileIn);

    // Creating the sequence to index
    DnaString genome;
    for (unsigned i = 0; i < length(ids); ++i){
        append(genome,seqs[i]);
    }
    // Parsing kmerFile
    StringSet<CharString> kmIds; // not used, but needed
    StringSet<CharString> kmSeqs;
    SeqFileIn kmSeqFileIn(toCString(kmerFile));
    readRecords(kmIds, kmSeqs, kmSeqFileIn);
    print("DONE");

    // Constants (shoulcound be in upper case)
    const int maxErr  = 2;               // Max number of errors, argument of the function
    const int nbRead  = length(ids);     // Number of reads in the file
    const int lenRead = length(seqs[0]); // size of the sequence to search in

    print("CREATING INDEX");
    Index<DnaString, BidirectionalIndex<FMIndex<TFastConfig> > > index(genome);
    print("DONE");
    // Result storage
    counter results;

    // Number of kmer, needed to conform with omp for loops
    int len = length(kmSeqs);
    // vector of Bool array, declared here so each thread have a private one
    std::array<std::vector<bool>, 3> tcount;
    omp_set_num_threads(nbThread);
    print("STARTING RESSEARCH");
    #pragma omp parallel firstprivate(nbRead,tcount,index)
    {
        auto delegateParallel = [&tcount, &nbRead, &lenRead](auto & iter, const DnaString & needle, int errors)
        {
            
            for (auto occ : getOccurrences(iter)){
                // Testing if the occurence is on one read only
                // (and not between two reads)
                if(occ%lenRead + length(needle) <= lenRead )
                {
                    tcount[errors][occ/lenRead] = true;
                   
                }
                std::cout << omp_get_thread_num() << " " << occ << " " << needle << " " << errors << std::endl;
            }
        };

        #pragma omp for schedule(dynamic)
        for(int k=0; k<len; k++)
        {
            std::string km = toCString(kmSeqs[k]);
            
            // setting temp counter and result value to false 
            for(int i =0; i<3; i++)
            {
                std::fill(tcount[i].begin(), tcount[i].end(), false);
                #pragma omp critical
                results[km][i] =  std::vector<bool>(nbRead,false);
            }

            // ressearch, filling tcount
            find<0, maxErr>(delegateParallel, index, kmSeqs[k], EditDistance() );

            // filling current element
            for(int error=0; error<3; error++)
            {
                for(int i = 0; i < nbRead ; i++ )
                {   
                    #pragma omp critical
                    results[km][error][i] = tcount[error][i];
                }
                #pragma omp critical
                std::cout << km << "\t" << error << std::endl;
                printVec(results[km][error]);
            }
        }
    }
    print("DONE");
    print("EXPORTING DATA");
    // printing results in file        
    printCounter(results, outputFile);
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
        "o", "outFile", "path to the output file",
        seqan::ArgParseArgument::STRING, "output file"));

    // Parse command line.
    seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);

    // If parsing was not successful then exit with code 1 if there were errors.
    // Otherwise, exit with code 0 (e.g. help was printed).
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res == seqan::ArgumentParser::PARSE_ERROR;

    // Extract option values and print them.
    std::string kmFile;     // default Kmer file
    std::string output = "out.txt";     // output file
    unsigned nbThread = 4;  // default number of thread
    unsigned nbErr = 2;     // default number of errors
    getOptionValue(kmFile, parser, "kf");
    getOptionValue(nbThread, parser, "nt");
    getOptionValue(output, parser, "o");
    
    print("PARSING ARGUMENTS");
    std::string text;
    getArgumentValue(text, parser, 0);
    std::ofstream outputFile;
    outputFile.open (output);
    outputFile << text << " " << kmFile << " " << nbThread << std::endl;
    print("DONE");
    findAdapt(text, kmFile, nbThread, outputFile);
    outputFile.close();
    print("END");

    return 0;
}