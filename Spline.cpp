/*
 * Spline.cpp
 *
 *  Created on: Jul 22, 2022
 *      Author: sergey
 */

#include "Spline.h"

namespace Skasp {

SplinePoints::SplinePoints()
{

}

SplinePoints::SplinePoints(CView1d _x, CView1d _y)
{
    Set(_x, _y);
}

SplinePoints::SplinePoints(CView1d _x, CView1d _y, CView1d _w)
{
    size_t n = std::min({_x.size(), _y.size(), _w.size()});
    points.resize(0);
    points.reserve(n);

    for (size_t i = 0 ; i < n; ++i)
    {
        if (i == 0 || _x[i] > _x[i - 1])
            points.emplace_back(_x[i], _y[i], _w[i]);
    }
}

void SplinePoints::Set(CView1d _x, CView1d _y)
{
    size_t n = std::min(_x.size(), _y.size());
    points.resize(0);
    points.reserve(n);

    for (size_t i = 0 ; i < n; ++i)
    {
        if (i == 0 || _x[i] > _x[i - 1])
            points.emplace_back(_x[i], _y[i], 1.0);
    }

    if (!points.empty())
    {
        points.front().w = 100.0;
        points.back().w = 100.0;
    }
}

Spline::Spline(size_t nInternalNodes, size_t spline_degree)
    : mK(spline_degree)
    , mKfact(factorial(mK))
    , mG(nInternalNodes)
    , mA(0.0)
    , mB(0.0)
    , mKnots(mG + 2*(mK + 1), 0.0)
    , mCoeff(mG + mK + 1)
    , mKnotsBuff(mG + 2, 0.0)
{

}

void Spline::GetKnots(Eigen::VectorXd& knots) const
{
    //knots.assign(mKnots.begin() + mK, mKnots.begin() + mG + mK + 2);
    std::copy(mKnots.begin() + mK, mKnots.begin() + mG + mK + 2, knots.begin());
}

const double* Spline::GetKnotsPtr() const
{
    return &mKnots[mK];
}

void Spline::SetKnots(const Eigen::VectorXd& h_knots)
{
    //
    std::copy(h_knots.begin(), h_knots.end(), mKnotsBuff.begin());

    //
    std::sort(mKnotsBuff.begin(), mKnotsBuff.end());

    //
    mA = mKnotsBuff[0];
    mB = mKnotsBuff[mG + 1];

    //
    std::vector<double>::iterator it = mKnots.begin();
    std::fill_n(it, mK + 1, mA);
    it = std::copy(++mKnotsBuff.begin(), --mKnotsBuff.end(), it + mK + 1);
    std::fill_n(it, mK + 1, mB);
}

void Spline::SetCoefficients(const Eigen::VectorXd& _coeff)
{
    //mCoeff.assign(_coeff.data(), _coeff.data() + _coeff.size());
    std::copy(_coeff.begin(), _coeff.end(), mCoeff.begin());
}

int Spline::GetLeftNodeIndex(double x, int min_id) const
{
    if (x < mA || x > mB)
        return -1;

    int l = min_id;
    while (l < mG + mK && (mKnots[l] > x || mKnots[l + 1] <= x))
        l += 1;
    return l;
}

double Spline::BSpline(double x, int deg, int id_left) const
{
    size_t id_right = id_left + deg + 1;

    if (x < mKnots[id_left] || x > mKnots[id_right])
        return 0;

    if (deg == 0)
        return !equal(x, mKnots[id_right]) ? 1.0 : 0.0;

    // if there are k + 1 coincident points on the left side
    if (mKnots[id_left + deg] < mKnots[id_left + deg + 1])
    {
        int j = 0;
        while (j < deg && equal(mKnots[id_left + j], mKnots[id_left + j + 1]))
            ++j;

        if (j == deg)
        {
            double alpha = (mKnots[id_left + deg + 1] - x)/(mKnots[id_left + deg + 1] - mKnots[id_left]);
            return std::pow(alpha, deg);
        }
    }

    // if there are k + 1 coincident points on the right side
    if (mKnots[id_left] < mKnots[id_left + 1])
    {
        int j = 1;
        while (j <= deg && equal(mKnots[id_left + j], mKnots[id_left + j + 1]))
            ++j;

        if (j == deg + 1)
        {
            double alpha = Alpha(x, id_left, deg + 1);
            return std::pow(alpha, deg);
        }
    }


    int l = GetLeftNodeIndex(x, id_left);
    SetBuffer(deg + 1);

    mBuff[id_left - l + deg] = 1.0;
    for (int j = 1; j < deg + 1; ++j)
        for (int i = l; i >= l - deg + j; --i)
            Propagate(mBuff, Alpha(x, i, 1 + deg - j), i - l + deg);

    return mBuff[deg];
}

// evaluate all B-splines of degree deg < k on interval at given point
bool Spline::BSplines(double x, int deg, std::vector<double>& splines) const
{
    if (deg > mK)
        return false;

    int l = GetLeftNodeIndex(x);
    splines[deg] = 1.0;

    for (int r = 1; r < deg + 1; ++r)
    {
        size_t v = l - r + 1;
        double w2 = Alpha(x, v + r, -r);
        splines[deg - r] = w2 * splines[deg - r + 1];
        for (int i = deg - r + 1; i < deg; ++i)
        {
            double w1 = w2;
            v++;
            w2 = Alpha(x, v + r, -r);
            splines[i] = (1.0 - w1) * splines[i] + w2 * splines[i + 1];
        }
        splines[deg] *= (1.0 - w2);
    }
    return true;
}

double Spline::BSplineDerivative(double x, int l, int i, int der_degree) const
{
    if (der_degree == 0)
        return BSpline(x, l, i);

    if (l == 0)
        return 0;

    double spline = 0;
    double c1 = mKnots[i + 1] - mKnots[i];
    double c2 = mKnots[i + 1 + 1] - mKnots[i + 1];
    if (allow_division(c1))
        spline += BSplineDerivative(x, l - 1, i, der_degree - 1) / c1;
    if (allow_division(c2))
        spline -= BSplineDerivative(x, l - 1, i + 1, der_degree - 1) / c2;
    return l*spline;
}

// get value of built spline at given point
double Spline::GetValue(double x) const
{
    if (x < mA || x > mB)
        return 0;

    int l = GetLeftNodeIndex(x);

    if (l < 0)
        return 0;

    SetBuffer(mK + 1);
    // De Boor Algorithm
    for (int i = 0; i < mK + 1; ++i)
        mBuff[i] = mCoeff[i + l - mK];

    for (int j = 1; j < mK + 1; ++j)
        for (int i = l; i >= l - mK + j; --i)
            Propagate(mBuff, Alpha(x, i, 1 + mK - j), i - l + mK);

    return mBuff[mK];
}

// get value of derivative of built spline at point x
double Spline::GetValueDerivative(double x, int der_degree) const
{
    if (der_degree == 1)
        return GetValue(x);

    if (der_degree > mK)
        return 0;

    int l = GetLeftNodeIndex(x);

    // if point is out of [a, b]
    if (l < 0)
    {
        if (der_degree > 1)
            return 0;
        return -1;  // todo: replace with derivatives outside of knots
    }

    double alpha = 1;
    double spline = 0;

    for (int i = 0; i < der_degree; ++i)
        alpha *= mK - i;

    SetBuffer(mK + 1);
    for (int i = 0; i < mK + 1; ++i)
        mBuff[i] = mCoeff[i + l - mK];

    for (int j = 1; j < der_degree; ++j)
        for (int i = j; i >= l - mK + j; --i)
        {
            size_t ndx = i - l + mK;
            mBuff[ndx] = (mBuff[ndx] - mBuff[ndx - 1]) / (mKnots[i + mK - j + 1] - mKnots[i]);
        }

    for (int i = der_degree; i < mK + 1; ++i)
        spline += mBuff[i] * BSpline(x, mK - der_degree, l + i - mK);

    return alpha * spline;
}

// get difference between leading derivatives of B-splines on interval [λ_{i}, λ_{i+k+1}] in λ_{q-} and λ_{q+}
double Spline::GetLeadDerivativeDifference(int i, int q) const
{
    if (i < q - mK - 1 || i > q)
        return 0;

    double num = static_cast<double>((2*(mK%2) - 1)*mKfact)*(mKnots[i + mK + 1] - mKnots[i]);
    double den = 1;
    for (int j = i; j < i + mK + 2; ++j)
    {
        if (j != q)
            den *= mKnots[q] - mKnots[j];
    }

    return num / den;
}

// get derivative of difference between leading derivatives
double Spline::GetLeadDerDiffDerKnot(int i, int q, int l) const
{
    if (l < i || l > i + mK + 1)
        return 0;
    return GetLeadDerDiffDerKnot(GetLeadDerivativeDifference(i, q), i, q, l);
}

// get derivative of difference between leading derivatives, if this difference is counted
double Spline::GetLeadDerDiffDerKnot(double lddk, int i, int q, int l) const
{
    int ndx = i + mK + 1;

    if (l < i || l > ndx)
        return 0;

    if (l != i && l != q && l != ndx)
        return lddk / (mKnots[q] - mKnots[l]);

    int c = (2 * (mK % 2) - 1) * mKfact;
    double product = c / lddk;
    double total_sum = 0;

    if (q != i && q != ndx)
    {
        if (l == i)
            return lddk / (mKnots[q] - mKnots[i]) / (mKnots[ndx] - mKnots[i]) * (mKnots[ndx] - mKnots[q]);
        if (l == i + mK + 1)
            return lddk / (mKnots[q] - mKnots[ndx]) / (mKnots[ndx] - mKnots[i]) * (mKnots[q] - mKnots[i]);
        if (q == l)
        {
            product *= mKnots[ndx] - mKnots[i];
            for (int j = i; j < i + mK + 2; ++j)
            {
                if (j != q)
                    total_sum += product / (mKnots[q] - mKnots[j]);
            }
            product *= product;
            return -c * (mKnots[ndx] - mKnots[i]) * total_sum / product;
        }
    }
    else
    {
        for (int j = i + 1; j < ndx; ++j)
            total_sum += product / (mKnots[q] - mKnots[j]);
        return -c * total_sum / SQ(product);
    }
    return 0;
}

// get value of derivative dS/dλ
double Spline::GetValueDerivativeKnot(double x, int knot_id) const
{
    if (x < mA)
    {
        if (knot_id == mK + 1)
        {
            double numerator = -mK * (mCoeff[1] - mCoeff[0]) * (x - mA);
            double denominator = SQ(mKnots[knot_id] - mA);
            return numerator / denominator;
        }
        return 0;
    }

   if (x > mB)
   {
       if (knot_id == mG + mK)
       {
            double numerator = -mK * (mCoeff[mG + mK] - mCoeff[mG + mK - 1]) * (x - mB);
            double denominator = SQ(mKnots[knot_id] - mB);
            return numerator / denominator;
       }
       return 0;
   }

   if (x <= mKnots[knot_id - mK] || x >= mKnots[knot_id + mK])
       return 0;

   int l = GetLeftNodeIndex(x);
   if (l < 0)
       return 0;

   if (l >= knot_id)
       ++l;

   SetBuffer(mK + 1);
   // De Boor algorithm
   for (int i = 0; i < mK + 1; ++i)
   {
       if (i < knot_id - l || i > knot_id - l + mK)
           mBuff[i] = 0;
       else
       {
           mBuff[i] = mCoeff[i + l - mK - 1] - mCoeff[i + l - mK];
           if (i + l + 1 <= knot_id)
               mBuff[i] /= mKnots[i + l + 1] - mKnots[i + l - mK];
           else if (i <= knot_id)
               mBuff[i] /= mKnots[i + l] - mKnots[i + l - mK];
           else
               mBuff[i] /= mKnots[i + l] - mKnots[i + l - mK - 1];
       }
   }

   for (int j = 1; j < mK + 1; ++j)
   {
       for (int i = l; i >= l - mK + j; --i)
       {
           double alpha = 0.0;
           if (i + mK - j + 1 <= knot_id)
               alpha = Alpha(x, i, mK - j + 1);
           else if (i <= knot_id)
               alpha = Alpha(x, i, mK - j);
           else
               alpha = Alpha(x, i - 1, mK - j + 1);

           Propagate(mBuff, alpha, i - l + mK);
       }
   }
   return mBuff[mK];
}

double Spline::Alpha(double x, int i, int k) const
{
    return (x - mKnots[i]) / (mKnots[i + k] - mKnots[i]);
}

void Spline::Propagate(std::vector<double>& N, double alpha, int i) const
{
    N[i] = alpha * N[i] + (1.0 - alpha) * N[i - 1];
}

void Spline::SetBuffer(size_t n) const
{
    if (mBuff.size() < n)
        mBuff.resize(n, 0.0);
    std::fill_n(mBuff.begin(), n, 0.0);
}