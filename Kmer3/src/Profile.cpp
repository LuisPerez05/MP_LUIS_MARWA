/*
 * Metodología de la Programación: Kmer3
 * Curso 2023/2024
 */

/** 
 * @file Profile.cpp
 * @author Luis Pérez Velasco <luispv05@correo.ugr.es>
 * @author Marwa Dris Azhir <marwadrisazhir@correo.ugr.es>
 * 
 * Created on 16 November 2023, 14:15
 */


#include "Profile.h"
#include <cmath>
#include <fstream>

using namespace std;

const string Profile::MAGIC_STRING_T="MP-KMER-T-1.0";

    Profile::Profile(){
        _profileId = "unknown";
        _size = 0;
    }

    Profile::Profile(int size){
        if((size < 0) || (size > DIM_VECTOR_KMER_FREQ)){
            throw std::out_of_range("Profile(int size): (size < 0) || (size > DIM_VECTOR_KMER_FREQ)");
        }
        _profileId = "unknown";
        _size = size;
    }
 
    const std::string& Profile::getProfileId() const{
        return _profileId;
    }

    void Profile::setProfileId(const std::string& id){
        _profileId = id;
    }
    
    const KmerFreq& Profile::at(int index) const{
        if((index < 0) || (index >= _size)){
            throw std::out_of_range("const KmerFreq& Profile::at(int index) const:"
                                    "(index < 0) || (index > _size)");
        }
        
        return _vectorKmerFreq[index];
    } 

    KmerFreq& Profile::at(int index){
         if((index < 0) || (index >= _size)){
            throw std::out_of_range("const KmerFreq& Profile::at(int index) const:"
                                    "(index < 0) || (index > _size)");
        }
        
        return _vectorKmerFreq[index];
    }

    int Profile::getSize() const{
        return _size;
    }
    
    int Profile::getCapacity() const{
        return DIM_VECTOR_KMER_FREQ;
    }
    
    double Profile::getDistance(const Profile& otherProfile) const{
        
        if((_size == 0) || (otherProfile.getSize() == 0)){
            throw std::invalid_argument("getDistance(const Profile& otherProfile) const: "
                                        "(_size == 0) || (otherProfile.getSize() == 0)");
        }
        
        
        double distance = 0;
        const int NOT_FOUND = otherProfile.getSize();
        
        for(int pos = 0; pos < _size; pos++){
            int rank = otherProfile.findKmer(_vectorKmerFreq[pos].getKmer());
            if(rank == -1) rank = NOT_FOUND;
            distance = distance + abs(pos - rank);
        }
        
        distance = distance / (_size * otherProfile.getSize());
        
        return distance;
    }
    
    int Profile::findKmer(const Kmer& kmer, int initialPos, int finalPos) const{
        
        bool found = false;
    
        while((initialPos <= finalPos)&&(!found)){
            if(_vectorKmerFreq[initialPos].getKmer().toString() == kmer.toString()){
                found = true;
            }else initialPos++;
        }

        if(!found) initialPos = -1;

        return initialPos;
    }

    int Profile::findKmer(const Kmer& kmer) const{

        return findKmer(kmer, 0, _size-1);
    }

    std::string Profile::toString() const{
         
        std::string sol = _profileId + "\n";
        
        sol += std::to_string(_size) + "\n";
    
        for(int n = 0; n < _size; n++){
            sol += _vectorKmerFreq[n].toString() + "\n";
        }
                
        return sol;
    }

    void Profile::sort(){
        for(int left = 1; left < _size; left++){
            KmerFreq to_insert = _vectorKmerFreq[left];
            int pos = left;

            while((pos > 0) && (greaterKmerFreq(to_insert, _vectorKmerFreq[pos-1]))){
                _vectorKmerFreq[pos] = _vectorKmerFreq[pos-1];
                pos--;
            }

            _vectorKmerFreq[pos] = to_insert;
        }
    }
    
    void Profile::save(const char fileName[]) const{
        
        ofstream output;
        output.open(fileName);
        
        if(!output){
            output.close();
            throw std::ios_base::failure("save(const char fileName[]) const: file cannot be opened");
        }
        else{
            output << MAGIC_STRING_T << endl;
            output << toString();

            if(!output){
                output.close();
                cerr<< "std::ios_base::failure(save(const char fileName[]) const: an error occurs while reading from the file)";
            }
            output.close();
        }
        
        
    } 

    void Profile::load(const char fileName[]){
        
        _size = 0;
        
        std::string magic_string;
        int nElements;
        int frequency;
        std::string nucleotides;
        
        ifstream input;
        input.open(fileName);
        
        
        
        if(!input){
            input.close();
            throw std::ios_base::failure("load(const char fileName[]): file cannot be opened");
        }
        
        input >> magic_string;
            
        if(magic_string != MAGIC_STRING_T){
            input.close();
            throw std::invalid_argument("load(const char fileName[]): magic_string != MAGIC_STRING_T");
        }
 
        input.ignore();
        getline(input, _profileId);
        input >> nElements;
        
        if((nElements >= DIM_VECTOR_KMER_FREQ) || (nElements < 0)){
            input.close();
            throw std::out_of_range("load(const char fileName[]): (nElements >= DIM_VECTOR_KMER_FREQ) || (nElements < 0)");
        }
            
        for(int n = 0; n < nElements; n++){
            
            input >> nucleotides;
            input >> frequency;
                
            KmerFreq aux;
            aux.setFrequency(frequency);
            aux.setKmer(Kmer(nucleotides));
                
            append(aux);
                
        }
            
        if(!input){
            input.close();
            throw std::ios_base::failure("load(const char fileName[]): an error occurs while reading from the file");
        }
            
            
        input.close();
    }
 
    void Profile::append(const KmerFreq& kmerFreq){
        
        bool same = false;
        int n = 0;
        while((n < _size) && !same){
            same = equal(kmerFreq, _vectorKmerFreq[n]);
            n++;
        }
        
        if(same){
            n--;
            _vectorKmerFreq[n].setFrequency(_vectorKmerFreq[n].getFrequency() + kmerFreq.getFrequency());
        }else{
            if(_size == DIM_VECTOR_KMER_FREQ){
                throw std::out_of_range("append(const KmerFreq& kmerFreq): _size == DIM_VECTOR_KMER_FREQ");
            }
            _vectorKmerFreq[_size] = kmerFreq;
            _size++;
        }
    }
    

    void Profile::normalize(const std::string& validNucleotides){
        
           

        // Normalize kmer i
        Kmer aux(_vectorKmerFreq[0].getKmer());
        aux.normalize(validNucleotides);
        _vectorKmerFreq[0].setKmer(aux);
        

        int n_non_repeated = 1;
        
        // Loop to traverse and normalize each one of the kmers in array
        // Loop to traverse the kmers in array from position 1 to position nElements-1
        for(int pos = 1; pos < _size; pos++){

              // Normalize kmer i
            Kmer aux(_vectorKmerFreq[pos].getKmer());
            aux.normalize(validNucleotides);
            _vectorKmerFreq[pos].setKmer(aux);

            
              // index = Position of array[i].getKmer() in the subarray that begins
              //         at position 0 and ends at position i-1
            int index = n_non_repeated - 1;
            int found;

            found = findKmer( _vectorKmerFreq[pos].getKmer(), 0, n_non_repeated - 1);
            
            if(found == -1){
                SwapElementsProfile(pos, n_non_repeated);
                n_non_repeated++;
            }else {
                _vectorKmerFreq[index].setFrequency(_vectorKmerFreq[index].getFrequency() + _vectorKmerFreq[pos].getFrequency());
            }


        }

        _size = n_non_repeated;
    }

    void Profile::deletePos(int pos){
        if((_size <= pos) || (pos < 0)){
            throw std::out_of_range("DeletePosArrayKmerFreq(KmerFreq array[], int& nElements, int pos): "
                                    "(nElements < pos) || (pos < 0)");
        }
        for(int element = pos; element < (_size - 1); element++){
            SwapElementsProfile(element, (element + 1));
        }
        _size--;
    }
        
    void Profile::zip (bool deleteMissing, int lowerBound){
        
        int pos = 0;
        while(pos < _size){

            bool no_delete = true;

            if(deleteMissing){
                std::string cadena = _vectorKmerFreq[pos].getKmer().toString();
                
                if(cadena.find(Kmer::MISSING_NUCLEOTIDE) != -1) no_delete = false;
                           
            }

            if(_vectorKmerFreq[pos].getFrequency() <= lowerBound){
                no_delete = false;
            }

            if(!no_delete){
                deletePos(pos);
            }else pos++;
        }
    }

    void Profile::join(const Profile& profile){
        int elements = profile.getSize();
        
        for(int n = 0; n < elements; n++){
            append(profile.at(n));
        }
    }
    
    void Profile::SwapElementsProfile(int first, int second){
        if((first > _size)||(second > _size)){
            throw std::out_of_range("SwapElementsArrayKmerFreq(int first, int second): "
                                    "(first > _size)||(second > _size)");
        }

        KmerFreq aux = _vectorKmerFreq[first];
        _vectorKmerFreq[first] = _vectorKmerFreq[second];
        _vectorKmerFreq[second] = aux;
    }
