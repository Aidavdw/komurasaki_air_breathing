#pragma once

#include "2dArray.h"

// Abstract class that represents a saved value, that can later be accessed throught the python interface.
class TwoDimensionalArrayRecord
{
public:
    TwoDimensionalArrayRecord() :
        src_(nullptr)
    {}
    explicit TwoDimensionalArrayRecord(const TwoDimensionalArray* src) :
        src_(src)
    {}

    std::vector<TwoDimensionalArray> records;
    void SaveRecord();

protected:
    const TwoDimensionalArray* src_;

};