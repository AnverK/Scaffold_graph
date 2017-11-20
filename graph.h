#ifndef GRAPH_H
#define GRAPH_H

#include "sam_reader.h"
#include <vector>
#include <unordered_map>
#include <algorithm>

class Graph
{
    friend class DotWriter;
public:

    Graph(SAM_Reader &sam_l, SAM_Reader &sam_r);

    Graph();

    ~Graph();

    struct ReadShortInfo                        //структура для hash-map
    {
        int contig, type, pos, len;
    };

    /*struct Edge
    {
        int v1, v2;
    };*/

    int SetInsertSize(string path);                             //необходимо указать путь к файлу .his, иначе insert_size не изменится

    int SetInsertSize(int n);

    int GetInsertSize();

    void SetTypeToReadAndPutReadInHash(int flag_left_sam);

    int CheckReadInHash(int flag_left_sam);

    int CheckReadInHashForInsertSize(int flag_left_sam, vector<int> &insert_size);

    void PutReadInHash(int flag_left_sam);

    void AddEdge(int type_1, int type_2, int contig_1, int contig_2);

    int SearchOptimalMinimalWeight(int flag, int num);

    void SetNewLibrary(SAM_Reader &sam_l, SAM_Reader &sam_r);

    void SetMinimalWeight(int w);

    void SaveMatrix(string path);

    void LoadMatrix(string path);

protected:
    unordered_map <string, ReadShortInfo> hash_read;
    int insert_size;
    SAM_Reader *sam_left, *sam_right;
    int cur_type;                                        //константа, показывающая как именно рид лежит на конце контига (принимает значения от 0 до 7)
    int min_weight;
    vector < vector <int> > adj_matr;

};

#endif // GRAPH_H
