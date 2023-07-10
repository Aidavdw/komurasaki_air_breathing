#include "debug/2darrayserialiser.h"
#include <fstream>
#include <iostream>
#include <sstream>

bool SaveStringToFile(const std::string& str, const std::string& filename)
{
    std::fstream file;
    file.open(filename, std::ios::out);
    if (!file)
        throw std::runtime_error("Error while outputting dump file");
    file << str;
    file.close();
    return true;
}

bool DumpTwoDimensionalArrayToFileGrid(const TwoDimensionalArray& src, const std::string& filename, const int significantDigits)
{
    std::stringstream ss;
    ss.precision(significantDigits);
    
    for (int x = 0; x < src.nX; x++)
    {
        for (int y = 0; y < src.nY; y++)
        {
            ss << src.GetAt(x,y);
            if (y != src.nY-1)
                ss << ", ";
        }
        if (x != src.nX-1)
            ss << "\n";
    }

    return SaveStringToFile(ss.str(), filename);
}

bool DumpPrettyWithGhostStrings(const TwoDimensionalArray& src, const std::string& filename,
    const int significantDigits)
{
    // Using stringstream to not have to re-allocate string every time.
    std::stringstream ss;
    ss.precision(significantDigits);

    // Also go over the ghost cells, so there is 2*nGhost more entries in each axis.
    for (int y = -src.nGhostCells; y < src.nY + src.nGhostCells; y++)
    {
        for (int x = -src.nGhostCells; x < src.nX + src.nGhostCells; x++)
        {
            // Vertical Spacing: If we're crossing the line between ghost -> non-ghost (bottom) or non-ghost -> ghost (top), skip a line to make it easier to see with the eye.
            if (x == -src.nGhostCells && ((y == 0) || (y == src.nY) ) )
                ss << "\n";

            // Horizontal spacing: dito, but now for the ghost cells on the sides
            if ((x == 0) || (x == src.nX))
                ss << "           , ";
            
            ss << src.GetAtWithGhost(x,y);
            if (x != src.nX + src.nGhostCells -1)
                ss << ", ";
        }
        ss << "\n";
    }

    return SaveStringToFile(ss.str(), filename);
}
