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
#include <cmath>
#include <fstream>
#include <string>

using namespace std;

/**
 * DEFAULT_VALID_NUCLEOTIDES is a c-string that contains the set of characters
 * that will be considered as valid nucleotides. 

 * The constructor of the class KmerCounter uses this c-string as a 
 * default parameter. It is possible to use a different c-string if that
 * constructor is used with a different c-string
 */
const char* const KmerCounter::DEFAULT_VALID_NUCLEOTIDES="ACGT";

KmerCounter::KmerCounter(int k, std::string validNucleotides) : _k(k), _validNucleotides(validNucleotides) {
    _allNucleotides = std::string(1, Kmer::MISSING_NUCLEOTIDE) + _validNucleotides;
    initFrequencies();
}

KmerCounter::KmerCounter(const KmerCounter& orig) : _k(orig._k), _validNucleotides(orig._validNucleotides), _allNucleotides(orig._allNucleotides) {
    int numRows = getNumRows();
    int numCols = getNumCols();
    _frequency = new int*[numRows];
    for (int i = 0; i < numRows; i++) {
        _frequency[i] = new int[numCols];
        for (int j = 0; j < numCols; j++) {
            _frequency[i][j] = orig._frequency[i][j];
        }
    }
}

KmerCounter::~KmerCounter() {
    int numRows = getNumRows();
    for (int i = 0; i < numRows; i++) {
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
    return pow(getNumNucleotides(), _k);
}

int KmerCounter::getNumberActiveKmers() {
    int count = 0;
    int numRows = getNumRows();
    int numCols = getNumCols();
    for (int i = 0; i < numRows; i++) {
        for (int j = 0; j < numCols; j++) {
            if (_frequency[i][j] > 0) {
                count++;
            }
        }
    }
    return count;
}

std::string KmerCounter::toString() {
    std::string str = _allNucleotides + " " + std::to_string(_k) + "\n";
    int numRows = getNumRows();
    int numCols = getNumCols();
    for (int i = 0; i < numRows; i++) {
        for (int j = 0; j < numCols; j++) {
            str += std::to_string(_frequency[i][j]) + " ";
        }
        str += "\n";
    }
    return str;
}

void KmerCounter::increaseFrequency(Kmer kmer, int frequency) {
    int row, column;
    getRowColumn(kmer, row, column);
    if (row!= -1 && column!= -1) {
        _frequency[row][column] += frequency;
    } else {
        throw std::invalid_argument("Invalid kmer");
    }
}

KmerCounter KmerCounter::operator=(KmerCounter orig) {
    if (this!= &orig) {
        int numRows = getNumRows();
        int numCols = getNumCols();
        for (int i = 0; i < numRows; i++) {
            delete[] _frequency[i];
        }
        delete[] _frequency;
        _k = orig._k;
        _validNucleotides = orig._validNucleotides;
        _allNucleotides = orig._allNucleotides;
        _frequency = new int*[numRows];
        for (int i = 0; i < numRows; i++) {
            _frequency[i] = new int[numCols];
            for (int j = 0; j < numCols; j++) {
                _frequency[i][j] = orig._frequency[i][j];
            }
        }
    }
    return *this;
}

KmerCounter KmerCounter::operator+=(KmerCounter kc) {
    if (_k!= kc._k || _validNucleotides!= kc._validNucleotides) {
        throw std::invalid_argument("Incompatible KmerCounters");
    }
    int numRows = getNumRows();
    int numCols = getNumCols();
    for (int i = 0; i < numRows; i++) {
        for (int j = 0; j < numCols; j++) {
            _frequency[i][j] += kc._frequency[i][j];
        }
    }
    return *this;
}

void KmerCounter::calculateFrequencies(char* fileName) {
    initFrequencies();
    std::ifstream file(fileName);
    if (!file) {
        throw std::ios_base::failure("Error opening file");
    }
    std::string line;
    while (std::getline(file, line)) {
        for (int i = 0; i <= line.size() - _k; i++) {
            Kmer kmer(line.substr(i, _k));
            increaseFrequency(kmer);
        }
    }
    file.close();
}

Profile KmerCounter::toProfile() {
    int kmerFrequencies;
    int numRows = getNumRows();
    int numCols = getNumCols();
    for (int i = 0; i < numRows; i++) {
        for (int j = 0; j < numCols; j++) {
            if (_frequency[i][j] > 0) {
                Kmer kmer = getKmer(i, j);
                kmerFrequencies = _frequency[i][j];
            }
        }
    }
    return kmerFrequencies;
}

int KmerCounter::getNumRows() {
    return pow(getNumNucleotides(), _k - _k / 2);
}

int KmerCounter::getNumCols() {
    return pow(getNumNucleotides(), _k / 2);
}

int KmerCounter::getIndex(const std::string& kmer) const {
    int index = 0;
    for (char c : kmer) {
        if (_allNucleotides.find(c) == std::string::npos) {
            return -1;
        }
        index = index * _allNucleotides.size() + _allNucleotides.find(c);
    }
    return index;
}

std::string KmerCounter::getInvertedIndex(int index, int nCharacters) const {
    std::string invertedIndex;
    for (int i = 0; i < nCharacters; i++) {
        invertedIndex += _allNucleotides[index % _allNucleotides.size()];
        index /= _allNucleotides.size();
    }
    return invertedIndex;
}

void KmerCounter::getRowColumn(const Kmer& kmer, int& row, int& column) const {
    std::string firstHalf = kmer.getSubstring(0, _k / 2);
    std::string secondHalf = kmer.getSubstring(_k / 2, _k);
    row = getIndex(firstHalf);
    column = getIndex(secondHalf);
}

Kmer KmerCounter::getKmer(int row, int column) {
    if (row < 0 || row >= getNumRows() || column < 0 || column >= getNumCols()) {
        throw std::invalid_argument("Invalid row or column");
    }
    std::string firstHalf = getInvertedIndex(row, _k / 2);
    std::string secondHalf = getInvertedIndex(column, _k - _k / 2);
    return Kmer(firstHalf + secondHalf);
}

void KmerCounter::initFrequencies() {
    int numRows = getNumRows();
    int numCols = getNumCols();
    _frequency = new int*[numRows];
    for (int i = 0; i < numRows; i++) {
        _frequency[i] = new int[numCols];
        for (int j = 0; j < numCols; j++) {
            _frequency[i][j] = 0;
        }
    }
}

int KmerCounter::operator()(int row, int column) const {
    return _frequency[row][column];
}

int KmerCounter::operator()(int row, int column) {
    return _frequency[row][column];
}