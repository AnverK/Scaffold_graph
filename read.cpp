#include "read.h"

string Read::GetQName()
{
    return read.qName;
}
string Read::GetRName()
{
    return read.rName;
}
string Read::GetCigar()
{
    return read.cigar;
}
string Read::GetRNext()
{
    return read.rNext;
}
string Read::GetSeq()
{
    return read.seq;
}
string Read::GetQual()
{
    return read.qual;
}
string Read::GetOther()
{
    return read.other;
}
int Read::GetFlag()
{
    return read.flag;
}
int Read::GetPos()
{
    return read.pos;
}
int Read::GetMapQ()
{
    return read.mapQ;
}
int Read::GetPNext()
{
    return read.pNext;
}
int Read::GetTLen()
{
    return read.TLen;
}

