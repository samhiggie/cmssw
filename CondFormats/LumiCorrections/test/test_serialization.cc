#include "CondFormats/Serialization/interface/Test.h"

#include "../src/headers.h"

int main()
{
    testSerialization<LumiCorrections>();
    //testSerialization<std::vector<LumiCorrections>>();
    //testSerialization<std::vector<LumiCorrections::aCorrection>>();
}
