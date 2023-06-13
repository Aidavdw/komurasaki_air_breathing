#include "..\..\public\pythoninterface\record.h"

void TwoDimensionalArrayRecord::SaveRecord()
{
    // First validate the pointer
    if (!src_)
        throw std::runtime_error("Could not verify 2d array record pointer");
    records.push_back(*src_);
}

py::array_t<double> TwoDimensionalArrayRecord::AsNumpyArray(const int atIdx) const
{
    if (atIdx < 0 || atIdx > records.size())
        throw std::out_of_range("Cannot access record, as it is out of range.");
    
    return TwoDimensionalArrayToNumpyArray(records.at(atIdx));
}
