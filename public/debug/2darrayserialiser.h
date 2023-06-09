#pragma once
#include "2dArray.h"

bool SaveStringToFile(const std::string& str, const std::string& filename);  // Opens a file and writes a string into it.
bool DumpTwoDimensionalArrayToFileGrid(const TwoDimensionalArray& src, const std::string& filename, const int significantDigits=3);    // Dumps it as a really simple csv. No ghost cells.
bool DumpPrettyWithGhostStrings(const TwoDimensionalArray& src, const std::string& filename, const int significantDigits=3);    // Makes a human readable output with ghost cells, for debug purposes.