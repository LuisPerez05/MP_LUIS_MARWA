/*
 * Metodología de la Programación: Kmer5
 * Curso 2023/2024
 */

/** 
 * @file KmerCounter.cpp
 * @author Silvia Acid Carrillo <acid@decsai.ugr.es>
 * @author Andrés Cano Utrera <acu@decsai.ugr.es>
 * @author Luis Castillo Vidal <L.Castillo@decsai.ugr.es>
 * @author Javier Martínez Baena <jbaena@ugr.es>
 * 
 * Created on 22 December 2023, 10:00
 */

#include "KmerCounter.h"

using namespace std;

/**
 * DEFAULT_VALID_NUCLEOTIDES is a c-string that contains the set of characters
 * that will be considered as valid nucleotides. 

 * The constructor of the class KmerCounter uses this c-string as a 
 * default parameter. It is possible to use a different c-string if that
 * constructor is used with a different c-string
 */
const char* const KmerCounter::DEFAULT_VALID_NUCLEOTIDES="ACGT";

std::string KmerCounter::toString() const{
    string outputString = _allNucleotides + " " + to_string(_k) + "\n";
    
    for(int row=0; row<this->getNumRows(); row++){
        for(int col=0; col<this->getNumCols(); col++){
            outputString += to_string((*this)(row,col)) + " ";
        }
        outputString += "\n";
    }
    
    return outputString;
}

int KmerCounter::getIndex(const std::string& kmer) const{
    int index = 0;
    int base = 1;

    for (size_t i = 0; i < kmer.size(); i++) {
        size_t pos = _allNucleotides.find(kmer[kmer.size()-i-1]);
        if (pos == string::npos)
            return -1;
        index += pos * base;
        base *= _allNucleotides.size();
    }
    return index;
}

string KmerCounter::getInvertedIndex(int index, int nCharacters) const {
    string result(nCharacters, Kmer::MISSING_NUCLEOTIDE);

    for (int i = result.size(); i > 0; i--) {
        result[i - 1] = _allNucleotides[index % _allNucleotides.size()];
        index = index / _allNucleotides.size();
    }
    return result;
}

KmerCounter::KmerCounter(int k, string validNucleotides) {
    _k = k;
    _validNucleotides = validNucleotides;
    _allNucleotides = Kmer::MISSING_NUCLEOTIDE + _validNucleotides;
    _frequency = new int*[_allNucleotides.size()];
    for (int i = 0; i < _allNucleotides.size(); i++) {
        _frequency[i] = new int[_allNucleotides.size()];
        for (int j = 0; j < _allNucleotides.size(); j++) {
            _frequency[i][j] = 0;
        }
    }
}

KmerCounter::KmerCounter(KmerCounter orig) {
    _k = orig._k;
    _validNucleotides = orig._validNucleotides;
    _allNucleotides = orig._allNucleotides;
    _frequency = new int*[_allNucleotides.size()];
    for (int i = 0; i < _allNucleotides.size(); i++) {
        _frequency[i] = new int[_allNucleotides.size()];
        for (int j = 0; j < _allNucleotides.size(); j++) {
            _frequency[i][j] = orig._frequency[i][j];
        }
    }
}

KmerCounter::~KmerCounter() {
    for (int i = 0; i < _allNucleotides.size(); i++) {
        delete[] _frequency[i];
    }
    delete[] _frequency;
}

int KmerCounter::getNumNucleotides() {
    return _allNucleotides.size();
}

int KmerCounter::getK() {
    return _k;
}

int KmerCounter::getNumKmers() {
    return pow(_allNucleotides.size(), _k);
}

int KmerCounter::getNumberActiveKmers() {
    int count = 0;
    for (int i = 0; i < _allNucleotides.size(); i++) {
        for (int j = 0; j < _allNucleotides.size(); j++) {
            if (_frequency[i][j] > 0) {
                count++;
            }
        }
    }
    return count;
}

void KmerCounter::increaseFrequency(Kmer kmer, int frequency) {
    if (kmer.getNumNucleotides() != _k) {
        throw invalid_argument("Invalid kmer length");
    }
    for (int i = 0; i < _validNucleotides.size(); i++) {
        if (kmer.getNucleotide(i) != _validNucleotides[i] && kmer.getNucleotide(i) != Kmer::MISSING_NUCLEOTIDE) {
            throw invalid_argument("Invalid nucleotide in kmer");
        }
    }
    int row = getIndex(kmer.getFirstHalf());
    int col = getIndex(kmer.getSecondHalf());
    if (row == -1 || col == -1) {
        throw invalid_argument("Invalid kmer");
    }
    _frequency[row][col] += frequency;
}

KmerCounter KmerCounter::operator=(KmerCounter orig) {
    if (this != &orig) {
        for (int i = 0; i < _allNucleotides.size(); i++) {
            delete[] _frequency[i];
        }
        delete[] _frequency;

        _k = orig._k;
        _validNucleotides = orig._validNucleotides;
        _allNucleotides = orig._allNucleotides;
        _frequency = new int*[_allNucleotides.size()];
        for (int i = 0; i < _allNucleotides.size(); i++) {
            _frequency[i] = new int[_allNucleotides.size()];
            for (int j = 0; j < _allNucleotides.size(); j++) {
                _frequency[i][j] = orig._frequency[i][j];
            }
        }
    }
    return *this;
}

int KmerCounter::getFrequency(Kmer kmer) {
    if (kmer.getNumNucleotides() != _k) {
        throw invalid_argument("Invalid kmer length");
    }
    for (int i = 0; i < _validNucleotides.size(); i++) {
        if (kmer.getNucleotide(i) != _validNucleotides[i] && kmer.getNucleotide(i) != Kmer::MISSING_NUCLEOTIDE) {
            throw invalid_argument("Invalid nucleotide in kmer");
        }
    }
    int row = getIndex(kmer.getFirstHalf());
    int col = getIndex(kmer.getSecondHalf());
    if (row == -1 || col == -1) {
        throw invalid_argument("Invalid kmer");
    }
    return _frequency[row][col];
}

vector<Kmer> KmerCounter::getKmers() {
    vector<Kmer> kmers;
    for (int i = 0; i < _allNucleotides.size(); i++) {
        for (int j = 0; j < _allNucleotides.size(); j++) {
            if (_frequency[i][j] > 0) {
                string kmerStr = "";
                kmerStr += _allNucleotides[i];
                kmerStr += _allNucleotides[j];
                kmers.push_back(Kmer(kmerStr));
            }
        }
    }
    return kmers;
}

int KmerCounter::getIndex(string nucleotide) {
    int index = -1;
    for (int i = 0; i < _allNucleotides.size(); i++) {
        if (nucleotide == _allNucleotides[i]) {
            index = i;
            break;
        }
    }
    return index;
}

bool KmerCounter::operator==(const KmerCounter& other) {
    if (_k != other._k || _validNucleotides != other._validNucleotides || _allNucleotides != other._allNucleotides) {
        return false;
    }
    for (int i = 0; i < _allNucleotides.size(); i++) {
        for (int j = 0; j < _allNucleotides.size(); j++) {
            if (_frequency[i][j] != other._frequency[i][j]) {
                return false;
            }
        }
    }
    return true;
}

bool KmerCounter::operator!=(const KmerCounter& other) {
    return !(*this == other);
}

KmerCounter& KmerCounter::operator+=(const KmerCounter& other) {
    if (_k != other._k || _validNucleotides != other._validNucleotides || _allNucleotides != other._allNucleotides) {
        throw invalid_argument("Invalid kmer counters for addition");
    }
    for (int i = 0; i < _allNucleotides.size(); i++) {
        for (int j = 0; j < _allNucleotides.size(); j++) {
            _frequency[i][j] += other._frequency[i][j];
        }
    }
    return *this;
}

KmerCounter KmerCounter::operator+(const KmerCounter& other) {
    KmerCounter result(*this);
    result += other;
    return result;
}

void KmerCounter::clear() {
    for (int i = 0; i < _allNucleotides.size(); i++) {
        for (int j = 0;j < _allNucleotides.size(); j++) {
            _frequency[i][j] = 0;
        }
    }
}

void KmerCounter::calculateFrequencies(char* fileName) {
    clear();
    ifstream file(fileName);
    if (!file.is_open()) {
        throw ios_base::failure("Could not open file");
    }
    string line;
    string nucleotide;
    Kmer kmer;
    while (getline(file, line)) {
        for (char c : line) {
            if (c == '\t' || c == '\n') {
                if (nucleotide.length() > 0) {
                    kmer = Kmer(nucleotide.substr(0, _k), _validNucleotides);
                    increaseFrequency(kmer, 1);
                }
                nucleotide = "";
            } else {
                nucleotide += toupper(c);
            }
        }
    }
    file.close();
}

Profile KmerCounter::toProfile() {
    Profile profile(_getNumberActiveKmers());
    vector<Kmer> kmers = getKmers();
    for (Kmer kmer : kmers) {
        profile.append(KmerFreq(kmer, getFrequency(kmer)));
    }
    return profile;
}

int KmerCounter::getNumRows() {
    return pow(_allNucleotides.size(), _k / 2);
}

int KmerCounter::getNumCols() {
    return pow(_allNucleotides.size(), _k / 2);
}

void KmerCounter::getRowColumn(Kmer kmer, int& row, int& column) {
    row = getIndex(kmer.getFirstHalf());
    column = getIndex(kmer.getSecondHalf());
}

Kmer KmerCounter::getKmer(int row, int column) {
    string firstHalf = getInvertedIndex(row, _k / 2);
    string secondHalf = getInvertedIndex(column, _k / 2);
    return Kmer(firstHalf + secondHalf);
}

void KmerCounter::initFrequencies() {
    for (int i = 0; i < _allNucleotides.size(); i++) {
        for (int j = 0; j < _allNucleotides.size(); j++) {
            _frequency[i][j] = 0;
        }
    }
}

int KmerCounter::operator()(int row, int column) {
    if (row < 0 || row >= getNumRows() || column < 0 || column >= getNumCols()) {
        throw invalid_argument("Invalid row or column");
    }
    return _frequency[row][column];
}

int& KmerCounter::operator()(int row, int column) {
    if (row < 0 || row >= getNumRows() || column < 0 || column >= getNumCols()) {
        throw invalid_argument("Invalid row or column");
    }
    return _frequency[row][column];
}
