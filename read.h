#ifndef READ_H
#define READ_H

#include <iostream>

using namespace std;

class Read
{
    friend class SAM_Reader;

public:
    string GetQName();
    string GetRName();
    string GetCigar();
    string GetRNext();
    string GetSeq();
    string GetQual();
    string GetOther();
    int GetFlag();
    int GetPos();
    int GetMapQ();
    int GetPNext();
    int GetTLen();

    struct ReadStr
    {
        string qName, rName, cigar, rNext, seq, qual, other;
        int flag, pos, mapQ, pNext, TLen;
    };

protected:
    ReadStr read;
};


#endif // READ_H
