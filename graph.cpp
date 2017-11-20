#include "graph.h"

Graph::Graph(SAM_Reader &sam_l, SAM_Reader &sam_r)
{
   sam_left = &sam_l;
   sam_right = &sam_r;
   insert_size = sam_l.min_len;
   min_weight = 0;

   for (unsigned int i = 0; i < adj_matr.size(); i++)
   {
       adj_matr[i].clear();
   }
   adj_matr.clear();

   adj_matr.resize(2*sam_l.num_nodes);
   for (unsigned int i = 0; i < adj_matr.size(); i++)
   {
       adj_matr[i].resize(adj_matr.size());
   }

}

Graph::Graph()
{
    min_weight = 0;
    for (unsigned int i = 0; i < adj_matr.size(); i++)
    {
        adj_matr[i].clear();
    }
    adj_matr.clear();
}

int Graph::SetInsertSize(string path)
{
    long long unsigned int len, sum=0, cur_sum=0;
    ifstream FIn;

    if(path.find(".his") != path.length() - 4 )
        return insert_size;

    // Начало чтения файла

    FIn.open(path, ios_base::in | ios_base::binary);
    FIn.seekg(0, ios_base::end);
    len = FIn.tellg()/4;

    vector <int> histogram;

    histogram.clear();
    histogram.resize( len );
    FIn.seekg(0, ios_base::beg);

    int i=0;
    while( !FIn.eof() )
    {
        FIn.read( (char*) &histogram[i], 4);
        i++;
        //cout << sum << endl;
        sum += histogram[i-1];
    }
    // Конец чтения файла


    i=0;
    while (cur_sum<sum/2)                //поиск медианы
    {
        cur_sum+=histogram[i];
        i++;
    }

    insert_size = i;

    return i;
}

int Graph::SetInsertSize(int n)
{
    insert_size = n;
    return n;
}

int Graph::GetInsertSize()
{
    return insert_size;
}

void Graph::SetTypeToReadAndPutReadInHash(int flag_left_sam)
{
    string q_name;
    int flag_forward_read = 0, flag_pos = -1;
    int cur_contig, pos;

    if(flag_left_sam)
    {
        q_name = sam_left->cur_read.GetQName();
        cur_contig = sam_left->GetCurContigNum();
        pos = sam_left->cur_read.GetPos();
        if(sam_left->cur_read.GetFlag() && 0x10)
            flag_forward_read = 1;
    }
    else
    {
        q_name = sam_right->cur_read.GetQName();
        cur_contig = sam_right->GetCurContigNum();
        pos = sam_right->cur_read.GetPos();
        if(sam_right->cur_read.GetFlag() && 0x10)
            flag_forward_read = 1;
    }


    if(q_name.length() <= 1)
        return;

    q_name.erase(q_name.end()-1);

    if(cur_contig < 0 )
    {
        hash_read.erase(q_name);
        return;
    }

    if(pos >= sam_left->GetContigLength(cur_contig) - insert_size)
        flag_pos = 1;
    else
        if(pos <= insert_size)
            flag_pos = 0;

    if(flag_pos == -1)
    {
        hash_read.erase(q_name);
        return;
    }

    cur_type = (flag_forward_read<<2) | (flag_pos<<1) | (flag_left_sam); 

    if( CheckReadInHash(flag_left_sam) == -1)
    {
        PutReadInHash(flag_left_sam);
        return;
    }
}

void Graph::AddEdge(int type_1, int type_2, int contig_1, int contig_2)
{
    if(contig_1 == contig_2)
        return;

    if(type_1 == 3 && type_2 == 4)
    {
        adj_matr[contig_1][contig_2]++;
        adj_matr[contig_2+sam_left->GetNumNodes()][contig_1+sam_left->GetNumNodes()]++;
        return;
    }

    if(type_1 == 4 && type_2 == 3)
    {
        adj_matr[contig_2][contig_1]++;
        adj_matr[contig_1+sam_left->GetNumNodes()][contig_2+sam_left->GetNumNodes()]++;
        return;
    }

    if(type_1 == 5 && type_2 == 4)
    {
        adj_matr[contig_1+sam_left->GetNumNodes()][contig_2]++;
        adj_matr[contig_2+sam_left->GetNumNodes()][contig_1]++;
        return;
    }

    if(type_1 == 4 && type_2 == 5)
    {
        adj_matr[contig_2+sam_left->GetNumNodes()][contig_1]++;
        adj_matr[contig_1+sam_left->GetNumNodes()][contig_2]++;
        return;
    }

    if(type_1 == 3 && type_2 == 2)
    {
        adj_matr[contig_1][contig_2+sam_left->GetNumNodes()]++;
        adj_matr[contig_2][contig_1+sam_left->GetNumNodes()]++;
        return;
    }

    if(type_1 == 2 && type_2 == 3)
    {
        adj_matr[contig_2][contig_1+sam_left->GetNumNodes()]++;
        adj_matr[contig_1][contig_2+sam_left->GetNumNodes()]++;
        return;
    }

    if(type_1 == 5 && type_2 == 2)
    {
        adj_matr[contig_2][contig_1]++;
        adj_matr[contig_1+sam_left->GetNumNodes()][contig_2+sam_left->GetNumNodes()]++;
        return;
    }

    if(type_1 == 2 && type_2 == 5)
    {
        adj_matr[contig_2][contig_1]++;
        adj_matr[contig_1+sam_left->GetNumNodes()][contig_2+sam_left->GetNumNodes()]++;
        return;
    }
}

int Graph::CheckReadInHash(int flag_left_sam)
{
    string q_name;
    int cur_contig;

    if(flag_left_sam)
    {
        q_name = sam_left->cur_read.GetQName();
        cur_contig = sam_left->GetCurContigNum();
    }
    else
    {
        q_name = sam_right->cur_read.GetQName();
        cur_contig = sam_right->GetCurContigNum();
    }

    q_name.erase(q_name.end()-1);

    unordered_map <string, ReadShortInfo> :: const_iterator got = hash_read.find(q_name);


    if(got == hash_read.end())
        return -1;

    AddEdge(cur_type, got->second.type, cur_contig, got->second.contig );

    hash_read.erase(q_name);
    return 1;
}

int Graph::CheckReadInHashForInsertSize(int flag_left_sam, vector <int> &insert_size)
{
    string q_name;
    int cur_contig;
    if(flag_left_sam)
    {
        q_name = sam_left->cur_read.GetQName();
        cur_contig = sam_left->GetCurContigNum();
    }
    else
    {
        q_name = sam_right->cur_read.GetQName();
        cur_contig = sam_right->GetCurContigNum();
    }

    if(q_name.size() > 1)
        q_name.erase(q_name.end()-1);
    else
        return 0;

    unordered_map <string, ReadShortInfo> :: const_iterator got = hash_read.find(q_name);

    if(got == hash_read.end())
        return 0;

    unsigned int size = 0;
    if(flag_left_sam)
    {
        if(got->second.contig == cur_contig)
        {
            if(got->second.pos >= sam_left->cur_read.GetPos())
            {
                size = abs(int(got->second.pos + got->second.len - sam_left->cur_read.GetPos()));
            }
            else
                size = abs(int(sam_left->cur_read.GetPos() + sam_left->cur_read.GetSeq().length() - got->second.pos));
            hash_read.erase (got);
        }
    }
    else
    {
        if(got->second.contig == cur_contig)
        {
            if(got->second.pos >= sam_right->cur_read.GetPos())
            {
                size = abs(int(got->second.pos + got->second.len - sam_right->cur_read.GetPos()));
            }
            else
                size = abs(int(sam_right->cur_read.GetPos() + sam_right->cur_read.GetSeq().length() - got->second.pos));

            hash_read.erase (got);
        }
    }
    if(size == 0)
        return 1;

    if(size > insert_size.size())
    {
       insert_size.resize(size+1);
    }
    insert_size[size]++;

    return 1;
}

void Graph::PutReadInHash(int flag_left_sam)
{
    ReadShortInfo info;
    string q_name;

    if(flag_left_sam)
    {
        q_name = sam_left->cur_read.GetQName();
        info.contig = sam_left->GetCurContigNum();
        info.pos = sam_left->cur_read.GetPos();
        info.len = sam_left->cur_read.GetSeq().length();
    }
    else
    {
        q_name = sam_right->cur_read.GetQName();
        info.contig = sam_right->GetCurContigNum();
        info.pos = sam_right->cur_read.GetPos();
        info.len = sam_right->cur_read.GetSeq().length();
    }
    if(q_name.size() > 1)
        q_name.erase(q_name.end()-1);
    else
        return;
    info.type = cur_type;
    hash_read.insert({q_name, info});
}

int Graph::SearchOptimalMinimalWeight(int flag, int num)
{
    vector <int> mas(adj_matr.size()*adj_matr.size());
    for(unsigned int i = 0; i<adj_matr.size(); i++)
    {
        for(unsigned int j = 0; j < adj_matr.size(); j++)
        {
            if(i != j)
            {
                if(2 * adj_matr[i][j] < adj_matr[j][i])                     //два ребра между двумя вершинами это не круто
                    mas[i*adj_matr.size() + j] = 0;
                else
                    mas[i*adj_matr.size() + j] = adj_matr[i][j];
            }
            else
                mas[i*adj_matr.size() + j] = 0;
        }
    }
    sort(mas.begin(), mas.end());

    if(flag == -1)
        min_weight = max(mas[mas.size() - adj_matr.size() + 2], 1);
    if(flag == 1)
        min_weight = max(mas[mas.size() - num], 1);
    if(flag == 2)
    {
        int num_of_edges = adj_matr.size()*(adj_matr.size() - 1) * num / 100;
        min_weight = max(mas[mas.size() - num_of_edges], 1);
    }

    return min_weight;
}

void Graph::SaveMatrix(string path)
{
    ofstream matr(path);
    if(matr.is_open() == 0)
        return;
    for(unsigned int i = 0; i < adj_matr.size(); i++)
    {
        for(unsigned int j = 0; j < adj_matr.size(); j++)
        {
            matr << adj_matr[i][j] << " ";
        }
        matr << endl;
    }
    matr.close();
}

void Graph::LoadMatrix(string path)
{
    ifstream matr(path);
    if(matr.is_open() == 0)
        return;
    for(unsigned int i = 0; i < adj_matr.size(); i++)
    {
        for(unsigned int j = 0; j < adj_matr.size(); j++)
        {
            matr >> adj_matr[i][j];
        }
    }
    matr.close();
}


void Graph::SetMinimalWeight(int w)
{
    min_weight = w;
}

Graph::~Graph()
{}

void Graph::SetNewLibrary(SAM_Reader &sam_l, SAM_Reader &sam_r)
{
    sam_left = &sam_l;
    sam_right = &sam_r;
    insert_size = sam_l.min_len;

    for (unsigned int i = 0; i < adj_matr.size(); i++)
    {
        adj_matr[i].clear();
    }
    adj_matr.clear();

    adj_matr.resize(2*sam_l.num_nodes);
    for (unsigned int i = 0; i < adj_matr.size(); i++)
    {
        adj_matr[i].resize(adj_matr.size());
    }
    hash_read.clear();
}

