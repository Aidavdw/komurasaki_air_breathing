#include "cell_values_container.h"

#include <cassert>

#define EPSILON 1E-7


void CellValues::Validate() const
{
    assert(density >  EPSILON);
    assert(p > EPSILON);
    assert(h > EPSILON);
}
