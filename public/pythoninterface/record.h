#pragma once
#include "conversion_functions.h"

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
    py::array_t<double> AsNumpyArray(const int idx) const;

protected:
    const TwoDimensionalArray* src_;

};