#include "..\..\public\pythoninterface\record.h"

void TwoDimensionalArrayRecord::SaveRecord()
{
    // First validate the pointer
    if (!src_)
        throw std::runtime_error("Could not verify 2d array record pointer");
    records.push_back(*src_);
}

py::array_t<double> TwoDimensionalArrayRecord::AsNumpyArray(const int idx) const
{
    if (idx < 0 || idx > records.size())
        throw std::out_of_range("Cannot access record, as it is out of range.");
    
    return TwoDimensionalArrayToNumpyArray(records.at(idx));
}
