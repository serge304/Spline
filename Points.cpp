/*
 * Points.cpp
 *
 *  Created on: Sep 28, 2022
 *      Author: sergey
 */

#include "Points.h"

namespace Skasp {

void Points2d::reserve(size_t n)
{
    args.reserve(n);
    vals.reserve(n);
}

void Points2d::resize(size_t n)
{
    args.resize(n);
    vals.resize(n);
}

void Points2d::push_back(double arg, double val)
{
    args.push_back(arg);
    vals.push_back(val);
}

void Points2d::push_back(const Point2d& point)
{
    push_back(point.x(), point.y());
}

void Points2d::push(size_t i, double arg, double val)
{
    args[i] = arg;
    vals[i] = val;
}

void Points3d::reserve(size_t n)
{
    args.reserve(n);
    vals1.reserve(n);
    vals2.reserve(n);
}

void Points3d::resize(size_t n)
{
    args.resize(n);
    vals1.resize(n);
    vals2.resize(n);
}

void Points3d::push_back(double arg, double val1, double val2)
{
    args.push_back(arg);
    vals1.push_back(val1);
    vals2.push_back(val2);
}

void Points3d::push_back(const Point3d& point)
{
    push_back(point.x(), point.y(), point.z());
}

} // namespace Skasp { namespace Base {
