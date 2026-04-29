/*
 * Math.cpp
 *
 *  Created on: Jun 20, 2023
 *      Author: sergey
 */

#include "Math.h"

namespace Skasp {

unsigned int factorial(unsigned int k)
{
    static const int table[] = {
            1, 1, 2, 6, 24, 120, 720, 5040, 40320, 362880, 3628800, 39916800, 479001600,
        };
    return table[k];
}

} // namespace Skasp {
