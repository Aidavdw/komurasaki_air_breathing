#include "..\..\public\pythoninterface\record.h"

void TwoDimensionalArrayRecord::SaveRecord()
{
    // First validate the pointer
    if (!src_)
        throw std::runtime_error("Could not verify 2d array record pointer");
    records.push_back(*src_);
}
