/*
 * Metodología de la Programación: Kmer5
 * Curso 2023/2024
 */

/** 
 * @file Profile.cpp
 * @author Silvia Acid Carrillo <acid@decsai.ugr.es>
 * @author Andrés Cano Utrera <acu@decsai.ugr.es>
 * @author Luis Castillo Vidal <L.Castillo@decsai.ugr.es>
 * @author Javier Martínez Baena <jbaena@ugr.es>
 * 
 * Created on 22 December 2023, 10:00
 */

#include "Profile.h"

using namespace std;

const string Profile::MAGIC_STRING_T="MP-KMER-T-1.0";
const string Profile::MAGIC_STRING_B="MP-KMER-B-1.0";

    Profile::Profile():_vectorKmerFreq (nullptr){
        _profileId = "unknown";
        _size = 0;
        allocate(INITIAL_CAPACITY);
    }

    Profile::Profile(int size):_vectorKmerFreq (nullptr){
        if(size < 0){
            throw std::out_of_range("Profile(int size): size < 0");
        }
        _profileId = "unknown";
        _size = size;
        allocate(size);
    }
    
    Profile::Profile(const Profile& orig):_vectorKmerFreq (nullptr){
        
        _profileId = orig.getProfileId();
        int nElements = orig.getSize();
        _size = nElements;
        allocate(orig.getCapacity());
        
        for(int pos = 0; pos < nElements; pos++){
            _vectorKmerFreq[pos] = orig.at(pos);
        }
    }
    
    Profile::~Profile(){
        deallocate();
    }
    
    Profile& Profile::operator=(const Profile& orig){
        if(&orig != this){
            
            delete[] this->_vectorKmerFreq;
                        
            copy(orig);
        }
        
        return *this;
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
            throw std::out_of_range(" KmerFreq& Profile::at(int index):"
                                    "(index < 0) || (index > _size)");
        }
        
        return _vectorKmerFreq[index];
    }

    int Profile::getSize() const{
        return _size;
    }
    
    int Profile::getCapacity() const{
        return _capacity;
    }
    
    double Profile::getDistance(const Profile& otherProfile) const{
        
        if((_size == 0) || (otherProfile.getSize() == 0)){
            throw std::invalid_argument("double Profile::getDistance(const Profile& otherProfile) const: "
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
            throw std::ios_base::failure("void Profile::save(const char fileName[]) const: file cannot be opened");
        }
        else{
            output << MAGIC_STRING_T << endl;
            output << toString();

            if(!output){
                output.close();
                cerr<< "std::ios_base::failure(void Profile::save(const char fileName[]) const: an error occurs while reading from the file)";
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
            throw std::ios_base::failure("void Profile::load(const char fileName[]): file cannot be opened");
        }
        
        input >> magic_string;
            
        if(magic_string != MAGIC_STRING_T){
            input.close();
            throw std::invalid_argument("void Profile::load(const char fileName[]): magic_string != MAGIC_STRING_T");
        }
 
        input.ignore();
        getline(input, _profileId);
        input >> nElements;
        
        if(nElements < 0){
            input.close();
            throw std::out_of_range("void Profile::load(const char fileName[]): nElements < 0");
        }
        
        allocate(nElements);
            
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
            throw std::ios_base::failure("void Profile::load(const char fileName[]): an error occurs while reading from the file");
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
            if(_size == _capacity){
                reallocate(_size + BLOCK_SIZE);
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
            throw std::out_of_range("void Profile::deletePos(int pos): "
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
            throw std::out_of_range("void Profile::SwapElementsArrayKmerFreq(int first, int second): "
                                    "(first > _size)||(second > _size)");
        }

        KmerFreq aux = _vectorKmerFreq[first];
        _vectorKmerFreq[first] = _vectorKmerFreq[second];
        _vectorKmerFreq[second] = aux;
    }
    
    void Profile::allocate(int capacity){
        delete[] _vectorKmerFreq;
        _vectorKmerFreq = new KmerFreq [capacity];
        _capacity = capacity;
    }
    
    void Profile::deallocate(){
        delete[] _vectorKmerFreq;
    }
    
    void Profile::reallocate(int newCapacity){
        KmerFreq* aux;
        aux = new KmerFreq [newCapacity];
        
        for(int pos = 0; pos < _size; pos++){
            aux[pos] = _vectorKmerFreq[pos];
        }
        
        delete[] _vectorKmerFreq;
        _vectorKmerFreq = aux;
        _capacity = newCapacity;
    }
    
    void Profile::copy(const Profile& orig){
        _size = orig.getSize();
        for(int pos = 0; pos < _size; pos++){
            _vectorKmerFreq[pos] = orig.at(pos);
        }
    }
    
    void Profile::write(std::ostream& outputStream) {
    outputStream << MAGIC_STRING_T << "\n";
    outputStream << _profileId << "\n";
    outputStream << _size << "\n";
    for (int i = 0; i < _size; ++i) {
        _vectorKmerFreq[i].write(outputStream);
    }
}

void Profile::read(std::istream& inputStream) {
    std::string magicString;
    std::getline(inputStream, magicString);
    if (magicString != MAGIC_STRING_T) {
        throw std::invalid_argument("Invalid magic string in Profile file");
    }
    std::getline(inputStream, _profileId);
    int size;
    inputStream >> size;
    if (size < 0) {
        throw std::out_of_range("Invalid size in Profile file");
    }
    allocate(size);
    for (int i = 0; i < size; ++i) {
        _vectorKmerFreq[i].read(inputStream);
    }
}

bool operator>(KmerFreq kmerFreq1, KmerFreq kmerFreq2) {
    if (kmerFreq1.getFrequency() > kmerFreq2.getFrequency()) {
        return true;
    } else if (kmerFreq1.getFrequency() == kmerFreq2.getFrequency()) {
        return kmerFreq1.getKmer() < kmerFreq2.getKmer();
    } else {
        return false;
    }
}

bool operator<(KmerFreq kmerFreq1, KmerFreq kmerFreq2) {
    return !(kmerFreq1 >= kmerFreq2);
}

bool operator==(KmerFreq kmerFreq1, KmerFreq kmerFreq2) {
    return kmerFreq1.getFrequency() == kmerFreq2.getFrequency() && kmerFreq1.getKmer() == kmerFreq2.getKmer();
}

bool operator!=(KmerFreq kmerFreq1, KmerFreq kmerFreq2) {
    return !(kmerFreq1 == kmerFreq2);
}

bool operator<=(KmerFreq kmerFreq1, KmerFreq kmerFreq2) {
    return kmerFreq1 < kmerFreq2 || kmerFreq1 == kmerFreq2;
}

bool operator>=(KmerFreq kmerFreq1, KmerFreq kmerFreq2) {
    return kmerFreq1 > kmerFreq2 || kmerFreq1 == kmerFreq2;
}

std::ostream& operator<<(std::ostream& os, const Profile& profile) {
    profile.write(os);
    return os;
}

std::istream& operator>>(std::istream& is, Profile& profile) {
    profile.read(is);
    return is;
}

