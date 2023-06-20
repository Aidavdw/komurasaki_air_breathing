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

void SetDomainInitialValueFromNumpyArrays(Domain& domain, const py::array_t<double>& rho, const py::array_t<double>& p, const py::array_t<double>& u, const py::array_t<double>& v)
{
    FillTwoDimensionalArrayFromNumpy(domain.p.currentTimeStep, p);
    FillTwoDimensionalArrayFromNumpy(domain.rho.currentTimeStep, rho);
    FillTwoDimensionalArrayFromNumpy(domain.u.currentTimeStep, u);
    FillTwoDimensionalArrayFromNumpy(domain.v.currentTimeStep, v);
}


void FillTwoDimensionalArrayFromNumpy(TwoDimensionalArray& out, const py::array_t<double>& pyArray)
{
    const py::buffer_info bufferInfo = pyArray.request();
    if (bufferInfo.ndim != 2)
        throw std::logic_error("Cannot transform numpy array to TwoDimensionalArray as it does not have 2 dimensions.");

    // check if they're the same size
    int rows = static_cast<int>(bufferInfo.shape[0]);
    int columns = static_cast<int>(bufferInfo.shape[1]);
    if (out.nX != columns || out.nY != rows)
        throw std::runtime_error("Cannot populate TwoDimensionalArray [" + std::to_string(out.nX) + ", " + std::to_string(out.nY) + "] from numpy array [" + std::to_string(columns) + ", " + std::to_string(rows) + "], as they are not the same shape");

    // Considers it as a flattened array itself, so get total amount of elements to iter over
    size_t nTot = 1;
    for (const auto r: bufferInfo.shape) {
        nTot *= static_cast<size_t>(r);
    }
    
    auto ptr = static_cast<double *>(bufferInfo.ptr);
    for (int xIdx = 0; xIdx < out.nX; xIdx++)
    {
        for (int yIdx = 0; yIdx < out.nY; yIdx++)
        {
            // Assign the value at the pointer, and move it forwards.
            out(xIdx, yIdx) = *ptr++;
        }
    }
}

