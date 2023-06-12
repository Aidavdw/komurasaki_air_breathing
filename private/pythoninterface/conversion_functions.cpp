#include "pythoninterface/conversion_functions.h"

py::array_t<double> TwoDimensionalArrayToNumpyArray(const TwoDimensionalArray& src)
{
    // Note that due to how indexing works in numpy, it is first y, then x. This is different from what is used in the rest of the code of this program.
    py::array_t<double, py::array::c_style> arr({src.nY, src.nX});
    auto ra = arr.mutable_unchecked();
    
    for (int yIndex = 0; yIndex < src.nY; yIndex++)
    {
        for (int xIndex = 0; xIndex < src.nX; xIndex++)
        {
            ra(yIndex, xIndex) = src.GetAt(xIndex, yIndex);
        }
    }

    return arr;
}
