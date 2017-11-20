#include "dotwriter.h"

DotWriter::DotWriter(Graph &g, int num_lib)
{
    graph = &g;
    path_to_dot = ".//output.dot";
    dot_file.open(path_to_dot);
    cur_lib = num_lib;
    colors.clear();
    colors.resize(6);
    colors[0] = "<blue>";
    colors[1] = "<red>";
    colors[2] = "<orange>";
    colors[3] = "<green>";
    colors[4] = "<black>";
    colors[5] = "<cyan>";
}

DotWriter::DotWriter(Graph &g, int num_lib, string path)
{
    graph = &g;    
    path_to_dot = path;
    dot_file.open(path_to_dot);
    cur_lib = num_lib;
    colors.clear();
    colors.resize(6);
    colors[0] = "<blue>";
    colors[1] = "<red>";
    colors[2] = "<orange>";
    colors[3] = "<green>";
    colors[4] = "<black>";
    colors[5] = "<cyan>";
}

void DotWriter::write_header()
{
    graph->sam_left->Rewind();

    string header_string;
    char buf[400];
    int length_of_contig;
    int flagBeg = 0;
    int num_nodes = 0;
    ContigName bufSt;
    string id="", cov="";

    dot_file << "digraph graph_picture {\nnode[fontname=<Courier> penwidth=<1.8> ]" << endl;

    for(int i = 0; flagBeg != 2; i++)
    {
        graph->sam_left->sam_file.getline(buf, sizeof(buf));
        header_string.append(buf);

        if( int( header_string.find("LN:") )!= -1)           //LN стоит в строках с объявлением контига
        {
            unsigned int pos = header_string.find("LN:")+3;
            length_of_contig = 0;
            while(pos < header_string.length())
            {
                //cout << s[i] << " ";
                length_of_contig *= 10;
                length_of_contig += header_string[pos] - '0';
                pos++;
            }
            if(length_of_contig >= graph->sam_left->min_len)
            {
                id.clear();
                cov.clear();
                num_nodes++;
                if( int ( header_string.find("ID_") ) != -1)
                {
                    int pos=header_string.find("ID_")+3;
                    while(header_string[pos]>='0' && header_string[pos]<='9')
                    {
                        id+=header_string[pos];
                        pos++;
                    }
                }
                if( int( header_string.find("cov_") ) != -1)
                {
                    int pos=header_string.find("cov")+4;
                    while( (header_string[pos]>='0' && header_string[pos] <='9') || header_string[pos] == '.')
                    {
                        cov+=header_string[pos];
                        pos++;
                    }
                }
                bufSt.cov = cov;
                bufSt.id = id;
                bufSt.len = length_of_contig;
                if(id[0] >= '0' && id[0] <= '9')
                    dot_file << "vertex_" << num_nodes-1 << "[label=\"ID: " << id << "\\n Len: " << length_of_contig << "\\n Cov: " << cov;
                else
                    dot_file << "vertex_" << num_nodes-1 << "[label=\"ID: " << "vertex_" << num_nodes-1 << "\\n Len: " << length_of_contig << "\\n Cov: " << cov;
                dot_file << "\", style=<filled> ,color=<black> ,fillcolor=<white> ]" << endl;
                if(id[0] >= '0' && id[0] <= '9')
                    dot_file << "INVERTED_vertex_" << num_nodes-1 << "[label=\"ID: " << id << "\\n Len: " << length_of_contig << "\\n Cov: " << cov;
                else
                    dot_file << "INVERTED_vertex_" << num_nodes-1 << "[label=\"ID: " << "INVERTED vertex_" << num_nodes-1 << "\\n Len: " << length_of_contig << "\\n Cov: " << cov;
                dot_file << "\", style=<filled> ,color=<black> ,fillcolor=<white> ]" << endl;
            }
            header_string.clear();
        }
        if( int( header_string.find("VN:") ) != -1)           //VN употребляется дважды, один раз в последней строке Header'a.
                                        //Когда строка будет с этим символом второй раз -- прекращается чтение Header'a.
        {
            flagBeg++;
            header_string.clear();
        }
    }


}

void DotWriter::write_graph()
{
    int num_of_edge = 0;

    unsigned int num_nodes = graph->adj_matr.size()/2;

    for (unsigned int i=0; i<graph->adj_matr.size(); i++)                             //удаление более лёгких рёбер между двумя вершинами и петель
    {
        for(unsigned int j=0; j<graph->adj_matr.size(); j++)
        {
            if(graph->adj_matr[i][j] > graph->adj_matr[j][i] || i==j)
            {
                graph->adj_matr[j][i] = 0;
            }
            else
            {
                graph->adj_matr[i][j] = 0;
            }
        }
    }

    for (unsigned int i=0; i<graph->adj_matr.size(); i++)
    {
        for(unsigned int j=0; j<graph->adj_matr.size(); j++)
        {
            if(graph->adj_matr[i][j] >= graph->min_weight)
            {
                if(i >= num_nodes )
                    dot_file << "INVERTED_";
                dot_file << "vertex_" << i % num_nodes << "->";
                if(j >= num_nodes )
                    dot_file << "INVERTED_";
                dot_file<< "vertex_" << j % num_nodes;
                dot_file << "[label=\"ID: " << num_of_edge << "\\n Weight: " << graph->adj_matr[i][j];
                dot_file << "\\n Lib#: " << cur_lib << "\", color=" << colors[cur_lib] << " ]" << endl;
                num_of_edge++;
            }
        }
    }
}

void DotWriter::write_end()
{
    dot_file << "}";
}

DotWriter::~DotWriter()
{
    dot_file.close();
}

void DotWriter::SetNumOfLibrary(int num)
{
    cur_lib = num;
}
