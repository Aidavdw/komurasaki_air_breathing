#pragma once
#include "conversion_functions.h"
#include "2dArray.h"

// Records are the way how this program handles returning data. Making a record of a property means you are interested in it, so its value will be stored for every time step. Represents a saved value, that can later be accessed through the python interface.
class TwoDimensionalArrayRecord
{
public:
    TwoDimensionalArrayRecord() :
        src_(nullptr)
    {}
    explicit TwoDimensionalArrayRecord(const TwoDimensionalArray* src) :
        src_(src)
    {}

    std::vector<TwoDimensionalArray> records; // All the recorded values for the source up to this point.
    void SaveRecord(); // Takes a snapshot of the current value of the source and adds it to the records array.
    py::array_t<double> AsNumpyArray(const int atIdx) const;  // returns the two dimensional array at a specific time step as a numpy array. Creates a new copy.

protected:
    const TwoDimensionalArray* src_;    // What this record is recording.

};