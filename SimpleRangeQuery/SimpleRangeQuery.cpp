// you can use includes, for example:
// #include <algorithm>

// you can write to stdout for debugging purposes, e.g.
// cout << "this is a debug message" << endl;
#include <algorithm>
#include <unordered_map>
#include <cmath>
#include <iostream>
#include <string>

const std::unordered_map<char, int> dnaValues{
    { 'A', 1 },
    { 'C', 2 },
    { 'G', 3 },
    { 'T', 4 },
};

std::vector<std::vector<int>> buildSparseTable(const std::string& seq) {

    int dnaSequenceLength = seq.size();
    int log2n = log2(dnaSequenceLength);
    std::vector<std::vector<int>> sparseTable(dnaSequenceLength, std::vector<int>(log2n + 1));
    for (uint16_t i = 0; i < dnaSequenceLength; ++i) {
        sparseTable[i][0] = i;
    }

    for (uint16_t j = 1; j <= log2n; j++) {
        for (int i = 0; i + (1 << j) <= dnaSequenceLength; ++i) {
            uint16_t left = sparseTable[i][j - 1];
            uint16_t right = sparseTable[i + (1 << (j - 1))][j - 1];
            sparseTable[i][j] = (dnaValues.at(seq[left]) <= dnaValues.at(seq[right])) ? left : right;
        }
    }

    return sparseTable;
}

int rangeMinimumQuery(const std::string& arr, int left, int right, const std::vector<std::vector<int>>& sparseTable) {
    int log2len = log2(right - left + 1);
    int left = sparseTable[left][log2len];
    int right = sparseTable[right - (1 << log2len) + 1][log2len];
    return std::min(dnaValues.at(arr[left]), dnaValues.at(arr[right]));
}

std::vector<int> processQueries(std::string& S, std::vector<int>& queriesStartPos, std::vector<int>& queriesEndPos) {
    std::vector<int> result(queriesStartPos.size(), 4);
    std::vector<std::vector<int>> table = buildSparseTable(S);

    for (uint16_t i = 0; i < queriesStartPos.size(); ++i) {
        result.at(i) = rangeMinimumQuery(S, queriesStartPos[i], queriesEndPos[i], table);
    }

    return result;
}
