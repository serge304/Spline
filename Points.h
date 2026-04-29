/*
 * Points.h
 *
 *  Created on: Sep 28, 2022
 *      Author: sergey
 */

#ifndef POINTS_H_
#define POINTS_H_

#include "Eigen.h"
#include <stddef.h>
#include <vector>

#if __cplusplus >= 202002L
#include <span>
#else
#include <span.h>
#endif

namespace Skasp {

//////////////////////////////////////////////////////////////////////////
///                                 Теги                               ///
//////////////////////////////////////////////////////////////////////////

enum
{
    ptArg = 0,
    ptVal = 1,
    ptVal1 = 1,
    ptVal2 = 2
};

//////////////////////////////////////////////////////////////////////////
///                               PointXd                              ///
//////////////////////////////////////////////////////////////////////////

typedef          double Point1d;
typedef Eigen::Vector2d Point2d;
typedef Eigen::Vector3d Point3d;

//////////////////////////////////////////////////////////////////////////
///                               PointsMap                            ///
//////////////////////////////////////////////////////////////////////////

typedef Eigen::Map<Eigen::VectorXd>        Map1d;
typedef Eigen::Map<const Eigen::VectorXd> CMap1d;

//////////////////////////////////////////////////////////////////////////
///                              PointsXd                              ///
//////////////////////////////////////////////////////////////////////////

typedef std::vector<Point1d>     Points1d;
typedef std::span<Point1d>       View1d;
typedef std::span<const Point1d> CView1d;

struct PointsBase
{
    //
    Points1d args;

    //
    PointsBase(size_t n = 0) : args(n) {}

    //
    size_t size() const { return args.size(); }
    bool  empty() const { return args.empty(); }

    //
    Map1d  arg_map()       { return  Map1d(&args[0], args.size()); }
    CMap1d arg_map() const { return CMap1d(&args[0], args.size()); }
};

struct Points2d : public PointsBase
{
    //
    Points1d vals;

    //
    Points2d(size_t n = 0) : PointsBase(n), vals(n) {}

    //
    void reserve(size_t n);
    void resize(size_t n);
    void push_back(double arg, double val);
    void push_back(const Point2d& point);
    void push(size_t i, double arg, double val);
};

struct Points3d : public PointsBase
{
    //
    Points1d vals1;
    Points1d vals2;

    //
    Points3d(size_t n = 0) : PointsBase(n), vals1(n), vals2(n) {}

    //
    void reserve(size_t n);
    void resize(size_t n);
    void push_back(double arg, double val1, double val2);
    void push_back(const Point3d& point);
};

template<typename T>
struct View2d_T
{
    std::span<T> args;
    std::span<T> vals;

    View2d_T(Points2d& points)
        : args(points.args)
        , vals(points.vals)
    {}

    View2d_T(const Points2d& points)
        : args(points.args)
        , vals(points.vals)
    {}

    View2d_T(Points2d& points, size_t n)
        : args(points.args.begin(), n)
        , vals(points.vals.begin(), n)
    {}

    View2d_T(const Points2d& points, size_t n)
        : args(points.args.begin(), n)
        , vals(points.vals.begin(), n)
    {}

    size_t size() const
    {
        return args.size();
    }
};

typedef View2d_T<double>        View2d;
typedef View2d_T<const double> CView2d;

} // namespace Skasp {

#endif /* POINTS_H_ */
