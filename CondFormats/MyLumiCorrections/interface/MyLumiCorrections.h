#ifndef MYLUMICORRECTIONS_H
#define MYLUMICORRECTIONS_H


/*
*Author: Sam Higginbotham
*Purpose: to save the corrections for the luminosity to the database 
*
*/
#include <sstream>
#include <cstring>
#include <vector>
#include <boost/serialization/vector.hpp>
#include "CondFormats/Serialization/interface/Serializable.h"

class MyLumiCorrections {

    public:
        //The types of corrections to save to the db
        //struct aCorrection { 
        //    float m_overallCorrection;
        //    float m_type1Fraction;
        //    float m_type1Residual;
        //    float m_type2Residual;
        //    std::vector<float> m_correctionsBX;
        //    COND_SERIALIZABLE;
        //};

        float m_overallCorrection;
        float m_type1Fraction;
        float m_type1Residual;
        float m_type2Residual;
        std::vector<float> m_correctionsBX;
        COND_SERIALIZABLE; 
        //std::vector<aCorrection> m_corrections;

};
#endif 
