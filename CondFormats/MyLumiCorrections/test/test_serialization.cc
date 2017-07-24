#include "CondFormats/Serialization/interface/Test.h"

#include "../src/headers.h"

int main()
{
    testSerialization<MyLumiCorrections>();
    //testSerialization<std::vector<MyLumiCorrections>>();
    //testSerialization<std::vector<MyLumiCorrections::aCorrection>>();
}
