/*
 * Spline.h
 *
 *  Created on: Jul 22, 2022
 *      Author: sergey
 */

#ifndef SPLINE_H_
#define SPLINE_H_

#include <Points.h>

namespace Skasp {

struct SplinePoint
{
    SplinePoint() : x(0.0), y(0.0), w(0.0) { }
    SplinePoint(double _x, double _y, double _w) : x(_x), y(_y), w(_w) { }
    double x, y, w;
};

class SplinePoints
{
public:
    SplinePoints();
    SplinePoints(CView1d _x, CView1d _y);
    SplinePoints(CView1d _x, CView1d _y, CView1d _w);

    size_t size() const { return points.size(); }
    bool empty() const { return points.empty(); }

    void Set(CView1d _x, CView1d _y);

    const SplinePoint& operator[](size_t ndx) const { return points[ndx]; }
    const SplinePoint& front() const { return points.front(); }
    const SplinePoint& back() const { return points.back(); }

public:
    std::vector<SplinePoint> points;
};

class Spline
{
public:
    Spline(size_t nInternalNodes, size_t spline_degree);

public:
    int GetDegree() const { return mK; }
    int GetInternalKnotsNum() const { return mG; }
    double GetLeftBound() const { return mA; }
    double GetRightBound() const { return mB; }
    int GetDegreeFactorial() const { return mKfact; }

    const double* GetCoefficients() const { return &mCoeff[0]; }
    void GetKnots(Eigen::VectorXd& knots) const;
    const double* GetKnotsPtr() const;

    void SetKnots(const Eigen::VectorXd& h_knots);
    void SetCoefficients(const std::vector<double>& _coeff) { mCoeff = _coeff; }
    void SetCoefficients(const Eigen::VectorXd& _coeff);
    void SetLeftEdge(double left_edge) { mA = left_edge; }
    void SetRightEdge(double right_edge) { mB = right_edge; }

public:
    //
    int GetLeftNodeIndex(double x, int min_id = 0) const;

    // get value of built spline at given point
    double GetValue(double x) const;

    // Evaluate all B-splines of degree deg < k on interval at given point
    bool BSplines(double x, int deg, std::vector<double>& splines) const;

    // get difference between leading derivatives of B-splines on interval [λ_{i}, λ_{i+k+1}] in λ_{q-} and λ_{q+}
    double GetLeadDerivativeDifference(int i, int q) const;

    // get value of derivative dS/dλ
    double GetValueDerivativeKnot(double x, int knot_id) const;

    // get derivative of difference between leading derivatives, if this difference is counted
    double GetLeadDerDiffDerKnot(double lddk, int i, int q, int l) const;

private:
    // Еvaluate B-spline of degree deg on interval [λ_{id}, λ_{id+deg+1}) at given point
    double BSpline(double x, int deg, int id) const;

    // Evaluate dD-th derivative of B-spline of degree l on interval [λ_i, λ_{i+l+1}) at given point
    double BSplineDerivative(double x, int l, int i, int der_degree = 1) const;

    // get value of derivative of built spline at point x
    double GetValueDerivative(double x, int der_degree = 1) const;

    // get derivative of difference between leading derivatives
    double GetLeadDerDiffDerKnot(int i, int q, int l) const;

private:
   double Alpha(double x, int i, int k) const;
   void Propagate(std::vector<double>& buff, double alpha, int i) const;
   void SetBuffer(size_t n) const;

private:
    int mK;                     // степень сплайна
    int mKfact;                 // факториал степени сплайна
    int mG;                     // количество внутренних точек
    double mA, mB;              // края диапазона для сплайна
    std::vector<double> mKnots; // вектор узлов сплайна
    std::vector<double> mCoeff; // вектор коэффициентов сплайна

    std::vector<double> mKnotsBuff;     // буфер для установки коэфициентов
    mutable std::vector<double> mBuff;  // буфер для расчетов
};

}  //namespace Skasp {


#endif /* SPLINE_H_ */
