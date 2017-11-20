#ifndef DotWriter_H
#define DotWriter_H

#include "graph.h"
#include <fstream>

class DotWriter
{
public:
    DotWriter(Graph &g, int num_lib);
    DotWriter(Graph &g, int num_lib, string path);
    ~DotWriter();
    void write_header();
    void write_graph();
    void SetNumOfLibrary(int num);

    struct ContigName
    {
        int len;
        string id, cov;
    };

    void write_end();

protected:
    Graph *graph;
    string path_to_dot;
    ofstream dot_file;
    int cur_lib;
    vector <string> colors;
};

#endif // DotWriter_H
