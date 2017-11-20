#ifndef SAM_READER_H_INCLUDED
#define SAM_READER_H_INCLUDED

#pragma once

#include <iostream>
#include <fstream>
#include <math.h>
#include <vector>
#include <unordered_map>

#include "read.h"

using namespace std;

class Read;

class SAM_Reader
{
    friend class Graph;
    friend class DotWriter;
public:
    ~SAM_Reader();
    void LoadLibrary(string path);
    void ReadHeader();
    int ReadNext();
    void Rewind();
    int GetNumNodes();
    int GetContigLength(int node_num);
    int GetContigNum(int node_num);
    int GetContigNum(string r_name);
    int GetContigNum(Read read);
    int GetCurContigNum();
    void SetContigMinLength(int num);
    ifstream sam_file;


protected:

    unordered_map <string, int> contig_num;      //пара из ReferenceName и номера контига в векторе из "хороших" контигов

    vector <int> contig_length;

    int flag_header;       // прочитан ли Header
    string path_to_sam_file;

    int num_nodes;
    Read cur_read;

    int min_len;

};
#endif // SAM_READER


