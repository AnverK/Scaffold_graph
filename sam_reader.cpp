#include "sam_reader.h"

void SAM_Reader::LoadLibrary(string path)
{
    if(path.find_last_of(".sam")==path.length()-1)
    {
        if(sam_file.is_open())
            sam_file.close();
        sam_file.open(path.c_str(),ios::in);
        sam_file.seekg(0,ios_base::end);
        sam_file.seekg(0,ios_base::beg);
        flag_header = 0;
        num_nodes = 0;
        min_len = 7000;
        path_to_sam_file = path;
    }
}

void SAM_Reader::ReadHeader()
{
    if(flag_header==1) return;
    if(!sam_file.is_open()) return;

    string header_string;

    char buf[400];

    int length_of_contig;

    int flagBeg = 0;

    num_nodes = 0;

    contig_length.clear();


    for(int i = 0; flagBeg != 2; i++)
    {

        sam_file.getline(buf, sizeof(buf));

        header_string.append(buf);


        if( int( header_string.find("LN:") ) != -1)           //LN стоит в строках с объявлением контига
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

            if(length_of_contig >= min_len)
            {
                contig_length.push_back(length_of_contig);
                num_nodes++;

                if( int( header_string.find("SN:") ) != -1)
                {
                   //cout << header_string << endl;
                    int beg = header_string.find("SN:")+3;
                    int end = header_string.find("\tLN:");
                    string hash_str;
                    hash_str = header_string.substr(beg,end-beg);
                    contig_num.insert({hash_str, num_nodes-1});
                //    cout << hash_str << endl;
                }

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

    flag_header = 1;
}

int SAM_Reader::ReadNext()
{
    if(flag_header==0) return 0;
    if(!sam_file.is_open()) return 0;

    int k=0;
    char charBuf='0';
    vector <string> splitted_string(12);

    charBuf = sam_file.get();

    while(k<12 && !sam_file.eof())
    {
        while(((charBuf!='\t' && charBuf!=' ') ||k==11) && charBuf!='\n' && !sam_file.eof())
        {
            splitted_string[k].push_back(charBuf);
            charBuf = sam_file.get();
        }
        if(!sam_file.eof())
            while(charBuf=='\t' || charBuf==' ')
                charBuf = sam_file.get();
        k++;
    }

    cur_read.read.qName = splitted_string[0];
    cur_read.read.flag = atoi(splitted_string[1].c_str());
    cur_read.read.rName = splitted_string[2];
    cur_read.read.pos = atoi(splitted_string[3].c_str());
    cur_read.read.mapQ = atoi(splitted_string[4].c_str());
    cur_read.read.cigar = splitted_string[5];
    cur_read.read.rNext = splitted_string[6];
    cur_read.read.pNext = atoi(splitted_string[7].c_str());
    cur_read.read.TLen = atoi(splitted_string[8].c_str());
    cur_read.read.seq = splitted_string[9];
    cur_read.read.qual = splitted_string[10];
    cur_read.read.other = splitted_string[11];


    if(sam_file.eof()) sam_file.close();


    return 1;
}

void SAM_Reader::Rewind()
{
    if(!sam_file.is_open())
        if(path_to_sam_file.find_last_of(".sam")==path_to_sam_file.length()-1)
        {
            sam_file.open(path_to_sam_file.c_str(),ios::in);
            flag_header=0;
            num_nodes=0;
        }
    sam_file.seekg(0,ios::beg);
}

int SAM_Reader::GetContigLength(int node_num)
{
    return contig_length.at(node_num);
}

int SAM_Reader::GetContigNum(string r_name)
{
    unordered_map <string, int> :: const_iterator got = contig_num.find(r_name);

    if(got == contig_num.end())
        return -1;

    return got->second;
}

int SAM_Reader::GetContigNum(Read read)
{
    return GetContigNum(read.GetRName());
}

int SAM_Reader::GetCurContigNum()
{
    return GetContigNum(cur_read.GetRName());
}

int SAM_Reader::GetNumNodes()
{
    return num_nodes;
}

SAM_Reader::~SAM_Reader()
{
    flag_header=0;
    if(sam_file.is_open())
    {
        Rewind();
        sam_file.close();
    }
    num_nodes=0;
}

void SAM_Reader::SetContigMinLength(int num)
{
    min_len = num;
}


