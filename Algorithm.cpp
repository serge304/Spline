/*
 * Algorithm.cpp
 *
 *  Created on: Jun 25, 2023
 *      Author: sergey
 */

#include "Algorithm.h"
#include <cmath>

namespace Skasp {

double round_one_dp(double x)
{
	int mul = 10;
	return x > 0.0 ?
			std::floor(x*mul + 0.5) / mul :
			std::ceil (x*mul - 0.5) / mul;
}

double round_two_dp(double x)
{
    int mul = 100;
    return x > 0.0 ?
               std::floor(x*mul + 0.05) / mul :
               std::ceil (x*mul - 0.05) / mul;
}

} // namespace Skasp {
