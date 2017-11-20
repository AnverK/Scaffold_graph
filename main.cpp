#include <iostream>
#include "sam_reader.h"
#include "graph.h"
#include "dotwriter.h"
#include "yaml-cpp/yaml.h"
#include <fstream>

using namespace std;


struct library
{
    string insert_size, left_sam, right_sam;
    int min_weight, num_of_edges, density_of_graph, flag, min_contig_len;
};

void count_IS(string path, SAM_Reader &s_l, SAM_Reader &s_r, Graph &g)
{
    vector <int> insert_size(1000);
    int k = 0, ok1 = 1, ok2 = 1;
    while( ok1 && ok2 )                         //читаем sam-файл до конца
    {
        ok1 = s_l.ReadNext();
        ok2 = s_r.ReadNext();
        if(g.CheckReadInHashForInsertSize(1, insert_size) == 0)
            g.PutReadInHash(1);

        if(g.CheckReadInHashForInsertSize(0, insert_size) == 0)
            g.PutReadInHash(0);
        if(k > 100000)
        {
            cout << s_l.sam_file.tellg()/1024/1024 << endl;
            k = 0;
        }
        k++;
    }
    while(ok1)
    {
        ok1 = s_l.ReadNext();
        if(g.CheckReadInHashForInsertSize(1, insert_size) == 0)
            g.PutReadInHash(1);
        if(k > 100000)
        {
            cout << s_l.sam_file.tellg()/1024/1024 << endl;
            k = 0;
        }
        k++;
    }
    while(ok2)
    {
        ok2 = s_r.ReadNext();
        if(g.CheckReadInHashForInsertSize(0, insert_size) == 0)
            g.PutReadInHash(0);
        if(k > 100000)
        {
            cout << s_l.sam_file.tellg()/1024/1024 << endl;
            k = 0;
        }
        k++;
    }
    ofstream FInsertSize;
    FInsertSize.open(path,ios_base::out|ios_base::binary|ios_base::trunc);
    FInsertSize.write((char *)&insert_size[0],insert_size.size()*4);
    FInsertSize.close();
    ok1 = 1;
    ok2 = 1;
    s_l.Rewind();
    s_r.Rewind();
    s_l.ReadHeader();
    s_r.ReadHeader();
    g.SetNewLibrary(s_l, s_r);
}

void read_library(SAM_Reader &s_l, SAM_Reader &s_r, Graph &g)
{
    int ok1 = 1, ok2 = 1, k = 0;
    while(ok1 && ok2)
    {
        ok1 = s_l.ReadNext();
        ok2 = s_r.ReadNext();
        g.SetTypeToReadAndPutReadInHash(1);
        g.SetTypeToReadAndPutReadInHash(0);
        k++;                                    // Просто показывает сколько Мб прочитано
        if(k == 100000)
        {
            cout << s_l.sam_file.tellg()/1024/1024 << endl;
            k=0;
            //break;
        }
    }
    while(ok1 )
    {
        ok1 = s_l.ReadNext();
        g.SetTypeToReadAndPutReadInHash(1);

    }

    while(ok2 )
    {
        ok2 = s_r.ReadNext();
        g.SetTypeToReadAndPutReadInHash(0);
    }
}



int main(int argc, char *argv[])
{
    YAML::Node config;
    if(argc == 2)
    {
        std::ifstream fin(argv[1]);
        YAML::Parser parser(fin);
        parser.GetNextDocument(config);
//        config
//        config = YAML::LoadFile(argv[1]);
    }
    else
    {
        cout << "You should write the path to YAML-file as a parameter" << endl;
        return 0;
    }

    string path_graph;
    int num_libs;
    config["output path"] >> path_graph;
    config["number of libraries"] >> num_libs;

    Graph g;
    SAM_Reader s_l, s_r;

    DotWriter dot(g, 0, path_graph);

    int flag_count_IS;

    library lib;

    for(int i = 0; i < num_libs; i++)
    {
        cout << "Preparing for making " << i+1 << " library..." << endl;

        lib.flag = -1;
        lib.min_weight = -1;
        lib.density_of_graph = -1;
        lib.num_of_edges = -1;
        lib.insert_size = "count";
        flag_count_IS = 0;
        lib.min_contig_len = 7000;

        dot.SetNumOfLibrary(i);

        int insert_size = -1;
        string buf;

        config["libraries"][i]["left sam"] >> lib.left_sam;
        config["libraries"][i]["right sam"] >> lib.right_sam;
        config["libraries"][i]["insert size"] >> buf;

        if(!buf.empty())
        {
            config["libraries"][i]["insert size"] >> lib.insert_size;
        }

        if(lib.insert_size == "count")
        {
            flag_count_IS = 1;
        }
        else
            if(lib.insert_size.find_last_of(".his") != lib.insert_size.length()-1)
            {
                int flag = 0;
                for(unsigned int j = 0; j < lib.insert_size.size(); j++)
                {
                    if(lib.insert_size[j] < '0' || lib.insert_size[j] > '9')
                    {
                        cout << "Warning: something has gone wrong with insert size. You could write only: \"count\" for "
                                "counting it (it takes as more time as reading library), number or path to .his-file" << endl;
                        cout << "Now programm would count it anyway" << endl;
                        flag = 1;
                        break;
                    }
                }
                if(flag == 0)
                     insert_size = atoi(lib.insert_size.c_str());
                else
                    flag_count_IS = 1;
            }

        buf.clear();
        config["libraries"][i]["minimal contig length"] >> buf;
        if(!buf.empty())
            config["minimal contig length"] >> lib.min_contig_len;

        s_l.LoadLibrary(lib.left_sam);
        cout << "left was written" << endl;
        s_l.SetContigMinLength(lib.min_contig_len);
        s_l.ReadHeader();
        cout << "header was written" << endl;


        s_r.LoadLibrary(lib.right_sam);
        cout << "right was written" << endl;
        s_r.SetContigMinLength(lib.min_contig_len);
        s_r.ReadHeader();
        cout << "header was written" << endl;

        g.SetNewLibrary(s_l, s_r);

        cout << "New Library was installed" << endl;

        if(flag_count_IS)
        {
            cout << "Counting Insert Size for " << i+1 << " library..." << endl;
            string buf = lib.left_sam;
            if( int(buf.find_last_of('\\')) != -1)
                buf.erase(buf.begin()+buf.find_last_of('\\')+1, buf.end());
            else
                buf.erase(buf.begin()+buf.find_last_of('/')+1, buf.end());
            lib.insert_size = buf + "insert_size_";
            lib.insert_size.push_back('0'+i);
            lib.insert_size.append(".his");
            count_IS(lib.insert_size, s_l, s_r, g);
            cout << "Counting Insert Size for " << i+1 << " library is over" << endl;

        }

        if(insert_size == -1)
            cout << "Insert Size: " << g.SetInsertSize(lib.insert_size) << endl;
        else
            cout << "Insert Size: " << g.SetInsertSize(insert_size) << endl;


        cout << endl << "Graph of " << i+1 << " library is making..." << endl;

        buf.clear();
        config["libraries"][i]["matrix"] >> buf;
        if(buf.empty())
        {
            read_library(s_l, s_r, g);
            string matrix_path = lib.left_sam;
            if( int(matrix_path.find_last_of('\\')) != -1)
                matrix_path.erase(matrix_path.begin()+matrix_path.find_last_of('\\')+1, matrix_path.end());
            else
                matrix_path.erase(matrix_path.begin()+matrix_path.find_last_of('/')+1, matrix_path.end());
            matrix_path = matrix_path + "matrix_";
            matrix_path.push_back('0'+i);
            matrix_path.append(".mtrx");
            g.SaveMatrix(matrix_path);
        }
        else {
            g.LoadMatrix(buf);
        }


        //Определяет количество рёбер в графе по одной из трёх характеристик (начало)

        buf.clear();
        config["libraries"][i]["minimal weight"] >> buf;
        if(!buf.empty())
        {
            config["libraries"][i]["minimal weight"] >> lib.min_weight;
            if(lib.min_weight > 0)
                lib.flag = 0;
        }
        buf.clear();
        config["libraries"][i]["number of edges"] >> buf;
        if(!buf.empty())
        {
            config["libraries"][i]["number of edges"] >> lib.num_of_edges;
            if(lib.num_of_edges > 0)
                lib.flag = 1;
        }
        buf.clear();
        config["libraries"][i]["density of graph"] >> buf;
        if(!buf.empty())
        {
            config["libraries"][i]["density of graph"] >> lib.density_of_graph;
            if(lib.density_of_graph > 0 && lib.density_of_graph <= 100)
                lib.flag = 2;
        }

        //Определяет количество рёбер в графе по одной из трёх характеристик (конец)

        cout << lib.flag << endl << endl;

        if(lib.flag == -1)
        {
            cout << "Minimal weight was counted automatically" << endl;
            g.SearchOptimalMinimalWeight(-1, 0);
        }
        if(lib.flag == 0)
        {
            cout << "Minimal weight was counted by... minimal weight" << endl;
            g.SetMinimalWeight(lib.min_weight);
        }
        if(lib.flag == 1)
        {
            cout << "Minimal weight was counted by... number of edges" << endl;
            g.SearchOptimalMinimalWeight(1, lib.num_of_edges);
        }
        if(lib.flag == 2)
        {
            cout << "Minimal weight was counted by... density of graph" << endl;
            g.SearchOptimalMinimalWeight(2, lib.density_of_graph);
        }



        cout << "Graph of " << i+1 << " library has been built." << endl;
        cout << "Dot file of " << i+1 << " library is making..." << endl;

       // DotWriter dot(g,0, "E:\\QTProjects\\Practice\\TESTS\\test1\\graph.dot");

        if(i == 0)
        {
            cout << "DOT Header" << endl;
            dot.write_header();             //ОСТОРОЖНО! в этой функции происходит Rewind sam_left
        }

        cout << "DOT Graph" << endl;
        dot.write_graph();


        cout << "Database of " << i+1 << " library was written in DOTfile." << endl;


        cout << i+1 << " library was maden" << endl;

    }
    cout << "the end" << endl;
    return 0;
}

