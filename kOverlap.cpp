#include<cstdlib>
#include<iostream>
#include<fstream>
#include <unordered_map>
#include <vector>
#include <array>
#include <algorithm>
#include <ctype.h>


using namespace std;

// Check if two strings have maximum overlap
string haveOverlap(string seq1, string seq2)
{   
    int l1   = seq1.size();
    int l2   = seq2.size();
    auto it1 = seq1.begin();
    auto it2 = seq2.begin();
    int minOverlap = min(l1,l2)-1 ;
    string overlap1 (it1+(l1-minOverlap), seq1.end());
    string overlap2 (it2, it2+minOverlap);
    
    if( overlap1 == overlap2)
    {
        string endSeq2 (it2+minOverlap,seq2.end());
        return(seq1 + endSeq2);
    }
    return("");
}

template <typename TContainer>
void splitLine(string text, string separator, TContainer &container)
{
    container.clear();
    string tempString;
    for(char c: text)
    {
        string s(1, c);
        if(s == separator)
        {
            container.push_back(tempString);
            tempString = "";
        }
        else if( isalnum(c))
        {
            tempString.append(s);
        }
    }
    container.push_back(tempString);
}

template <typename TKmCounter, typename TKmList>
void parseFile(string fileName ,string sep, TKmCounter &kmCount, TKmList &kmerList ){

    // input file stream
    ifstream inFile (fileName);

    if (inFile.is_open())
    {
        string line;
        string noerr, oneerr, twoerr, kmer;
        vector<string> lineSplit;

        getline (inFile,line); // skipping first line
        while ( getline (inFile,line) )
        {
            splitLine(line, sep, lineSplit);
            noerr  = lineSplit[0];
            oneerr = lineSplit[1];
            twoerr = lineSplit[2];
            kmer   = lineSplit[3];
            kmerList.push_back(kmer);
            kmCount[kmer] = {noerr, oneerr, twoerr};
        }
        inFile.close();
    }
    else{
        cout << "ERROR, COULD NOT OPEN FILE" << endl;
    }
    
}

template <typename TSmallVector>
string formatOutLine(string prefix, string key, TSmallVector errorCounts){
    string result = prefix + "\t" + key;
    for(auto elem: errorCounts){
        result += "\t" + elem;
    }
    return(result);
}


string basename(string path)
{   
    string currentName;
    for(char c : path){
        if(c == '/'){
            currentName = "";
        }
        else if(c!='.'){
            currentName += c;
        }
        else{
            return(currentName);        
        }
    }
    return("");
}

string  dirname(string path)
{
    string dir;
    string currentName;
    for(char c : path){
        if(c == '/'){
            dir += currentName + "/";
            currentName = "";
        }
        else if(c!='.'){
            currentName +=c;
        }
        else{
            return(dir);        
        }
    }
    return("");
}

int main(int argc, char** argv){

    // Initialising wariables
    // containers 
    vector<string> kmerList;
    unordered_map<string, array<string,3>> kmCount;
    // separator
    string sep = "\t";
    
    // Input file name 
    // string inputFile = argc>1 ? argv[1] : "example.txt";
    string inputFile = argv[1];
    //Creating outputfile using input file
    int l = inputFile.length();
    string outputFile = dirname(inputFile) + "overlap_output_"+ basename(inputFile) +".csv";
    parseFile(inputFile, sep, kmCount, kmerList);


    // searching longest overlap starting with the most frequent k-mer
    string km1 = kmerList[0];
    string ov = km1;
    // store element for final display
    vector<string> displayVector;
    // loop flag. Searching overlap untill we run out of k-mers or no overlap are found.
    bool found = true;
    displayVector.push_back( formatOutLine("FIRST", km1, kmCount[km1]) );    
    
    //     Avoid using the same k-mer twice.
    // /!\ Can lead to unpredictable behaviours with small k-mers
    vector<string> used = {km1};

    while(found){
        found = false;
        for(string km2: kmerList){
            if( std::find(used.begin(), used.end(), km2) == used.end() ){
                // searching right and left overlap
                string direct  = haveOverlap(ov, km2);
                string reverse = haveOverlap(km2, ov);
                if(direct != "" and reverse == ""){
                    ov = direct;
                    found = true;
                    used.push_back(km2);
                    displayVector.push_back(formatOutLine("RIGHT", km2 , kmCount[km2]));
                    break;
                }
                else if(reverse != "" and direct == ""){
                    ov = reverse;
                    found = true;
                    used.push_back(km2);
                    displayVector.insert(displayVector.begin(),formatOutLine("LEFT",  km2 , kmCount[km2]));
                    break;
                }

            }
        }
    }
    // printing in output file
    ofstream output (outputFile);
    output << ov << endl;
    for(auto line: displayVector){
        output << line << endl;
    }
    output.close();


    return(0);
}