/*
 * SplineController.cpp
 *
 *  Created on: Jul 25, 2022
 *      Author: sergey
 */

#include "SplineController.h"

namespace Skasp {

SplineController::SplineController()
  : mDelta(0.0)
  , mP(0.0)
  , mError(0.0)
  , mPenalty(0.0)
{

}

void SplineController::Reset()
{
    mDelta = 0.0;
    mP = 0.0;
    mError = 0.0;
    mPenalty = 0.0;

    mSplines.resize(0);
    mKnots.resize(0);
}

void SplineController::SetBuffers(int nInternalKnots, int degree)
{
    mKnots.resize(nInternalKnots + 2);
    mSplines.resize(degree + 1, 0.0);

    size_t dim = nInternalKnots + degree + 1;
    mA = Eigen::MatrixXd::Zero(dim, dim);
    mR = Eigen::VectorXd::Zero(dim);
    mC = Eigen::VectorXd::Zero(dim);
}

bool SplineController::Init(Spline& spline, const SplinePoints& points, int nInternalKnots)
{
    Reset();
    SetBuffers(nInternalKnots, spline.GetDegree());

    if (!InitiateGrid(spline, points, nInternalKnots))
        return false;

    Penalty(spline);
    return true;
}

bool SplineController::Init(Spline& spline, const SplinePoints& points)
{
    return Init(spline, points, spline.GetInternalKnotsNum());
}

double SplineController::Penalty(const Spline& spline)
{
    const double* knots = spline.GetKnotsPtr();
    int g = spline.GetInternalKnotsNum();
    mPenalty = 0.0;
    for (int i = 0; i < g + 1; ++i)
        mPenalty += 1.0 / (knots[i + 1] - knots[i]);
    mError = mDelta + mP * mPenalty;
    return mPenalty;
}

double SplineController::Delta(const Spline& spline, const SplinePoints& points, double smoothing_weight)
{
    mDelta = 0;
    for (size_t i = 0; i < points.size(); ++i)
    {
        const SplinePoint& point = points[i];
        double e = point.w * (point.y - spline.GetValue(point.x));
        mDelta += e * e;
    }

    double nu = 0;
    if (smoothing_weight > 0)
    {
        int k = spline.GetDegree();
        int g = spline.GetInternalKnotsNum();
        const double* coefficients = spline.GetCoefficients();
        for (int q = k + 1; q < g + k + 1; ++q)
        {
            double e = 0;
            for (int i = q - k - 1; i < q + 1; ++i)
            {
                e += coefficients[i] * spline.GetLeadDerivativeDifference(i, q);
                nu += e * e;
            }
        }
        nu *= smoothing_weight;
    }
    mDelta += nu;
    mError = mDelta + mP * mPenalty;
    return mDelta;
}

double SplineController::PenaltyDerivative(const Spline& spline, int knot_id)
{
    const double* knots = spline.GetKnotsPtr();
    double a = knots[knot_id - 1];
    double b = knots[knot_id];
    double c = knots[knot_id + 1];
    double bma = b - a;
    double cmb = c - b;
    return 1.0 / (cmb * cmb) - 1.0 / (bma * bma);
}

double SplineController::Error(const Spline& spline, const SplinePoints& points, double sw)
{
    return Delta(spline, points, sw) + mP * Penalty(spline);
}

double SplineController::ErrorGradient(const Spline& spline, const SplinePoints& points, double smoothing_weight, int knot_id)
{
    // least-square error
    double grad_error = 0;
    for (size_t i = 0; i < points.size(); ++i)
    {
        const SplinePoint& point = points[i];
        double w_sq = SQ(point.w);
        double diff = point.y - spline.GetValue(point.x);
        grad_error -= spline.GetValueDerivativeKnot(point.x, knot_id + spline.GetDegree()) * w_sq * diff;
    }

    // penalty error
    if (mP > 0)
        grad_error += 0.5 * mP * PenaltyDerivative(spline, knot_id);

    // smoothing error
    if (smoothing_weight > 0)
    {
        double sm_error = 0;
        int k = spline.GetDegree();
        int g = spline.GetInternalKnotsNum();
        const double* coefficients = spline.GetCoefficients();
        for (int q = k + 1; q < g + k + 1; ++q)
        {
            double sum1 = 0;
            double sum2 = 0;
            for (int i = q - k - 1; i < q + 1; ++i)
            {
                double ci = coefficients[i];
                double lead_der_diff = spline.GetLeadDerivativeDifference(i, q);
                sum1 += ci * lead_der_diff;
                sum2 += ci * spline.GetLeadDerDiffDerKnot(lead_der_diff, i, q, knot_id + spline.GetDegree());
            }
            sm_error += sum1 * sum2;
        }
        grad_error += smoothing_weight * sm_error;
    }
    return 2 * grad_error;
}

double SplineController::Theta(Spline& spline, const SplinePoints& points, double sw, double alpha, const Eigen::VectorXd& direction, const Eigen::VectorXd& fixed_knots)
{
    int g = spline.GetInternalKnotsNum();
    mKnots[0] = spline.GetLeftBound();
    mKnots[g + 1] = spline.GetRightBound();
    for (int i = 0; i < g; ++i)
         mKnots[i + 1] = fixed_knots[i + 1] + alpha * direction[i];
    spline.SetKnots(mKnots);

    if (Approximate(spline, points, sw))
        return Error(spline, points, sw);
    return -1;
}

bool SplineController::Approximate(Spline& spline, const SplinePoints& points, double smoothing_weight)
{
    int k = spline.GetDegree();
    int g = spline.GetInternalKnotsNum();
    size_t n = points.size();
    mA.setZero();
    mR.setZero();

    // spline error
    int l = 0;
    for (size_t r = 0; r < n; ++r)
    {
        const SplinePoint& point = points[r];
        l = spline.GetLeftNodeIndex(point.x, l);
        if (l < 0)
            return false;

        std::fill(mSplines.begin(), mSplines.end(), 0.0);
        if (!spline.BSplines(point.x, k, mSplines))
            return false;

        for (int i = 0; i < k + 1; ++i)
        {
            double w_sq = SQ(point.w);
            for (int j = 0; j < i + 1; ++j)
                mA(i + l - k, j + l - k) += w_sq * mSplines[i] * mSplines[j];
            mR(i + l - k) += w_sq * point.y * mSplines[i];
        }
    }

    // smoothing error
    if (smoothing_weight > 0)
    {
        for (int q = 0; q < g; ++q)
        {
            int ndx = q + k + 1;
            for (int i = q; i < ndx + 1; ++i)
            {
                double ai = spline.GetLeadDerivativeDifference(i, ndx);
                for (int j = q; j < i + 1; ++j)
                    mA(i,j) += smoothing_weight * ai * spline.GetLeadDerivativeDifference(j, ndx);
            }
        }
    }

    mC.noalias() = mA.selfadjointView<Eigen::Lower>().ldlt().solve(mR);
    spline.SetCoefficients(mC);

    return true;
}

bool SplineController::SpecDimensionalMinimization(Spline& spline, const SplinePoints& points, double sw,
                                   const Eigen::VectorXd& direction, const Eigen::VectorXd& error_derivative,
                                   const Eigen::VectorXd& fixed_knots)
{
    int g = spline.GetInternalKnotsNum();
    const double* knots = spline.GetKnotsPtr();
    double alpha_max = std::numeric_limits<double>::max(); //math.inf
    double a = spline.GetLeftBound();
    double b = spline.GetRightBound();

    if (direction[0] < 0)
        alpha_max = (a - knots[1]) / direction[0];

    for (int i = 0; i < g - 1; ++i)
    {
        if (direction[i] > direction[i + 1])
            alpha_max = std::min(alpha_max, (knots[i + 2] - knots[i + 1]) / (direction[i] - direction[i + 1]));
    }

    if (direction[g - 1] > 0)
        alpha_max = std::min(alpha_max, (b - knots[g]) / direction[g - 1]);

    double theta0 = mError;
    double theta0_der = direction.dot(error_derivative);
    double alpha0 = 0;
    double alpha2 = alpha_max / (1 - theta0 / alpha_max / theta0_der);
    double alpha1 = 0.5 * alpha2;
    double q0 = mDelta;
    double r0 = mPenalty;
    double theta1 = Theta(spline, points, sw, alpha1, direction, fixed_knots);
    if (theta1 < 0)
        return false;
    double q1 = mDelta;
    double r1 = mPenalty;

    int iteration = 0;
    int max_num_of_iterations = 1000;

    while (theta1 - theta0 >= 1.0e-4 && iteration < max_num_of_iterations)
    {
        double alpha_tilde = -0.5 * theta0_der * alpha1 * alpha1 / (theta1 - theta0 - theta0_der * alpha1);
        alpha1 = std::max(0.1 * alpha1, alpha_tilde);
        theta1 = Theta(spline, points, sw, alpha1, direction, fixed_knots);
        if (theta1 < 0)
            return false;
        q1 = mDelta;
        r1 = mPenalty;
        ++iteration;
    }

    if (iteration > 0)
    {
        if (theta1 > theta0)
            Theta(spline, points, sw, alpha0, direction, fixed_knots);
        // should we return false in the case of if?
        return true;
    }

    double theta2 = Theta(spline, points, sw, alpha2, direction, fixed_knots);
    if (theta2 < 0)
        return false;

    double q2 = mDelta;
    double r2 = mPenalty;

    while (theta2 < theta1)
    {
        alpha0 = alpha1;
        q0 = q1;
        r0 = r1;
        alpha1 = alpha2;
        theta1 = theta2;
        q1 = q2;
        r1 = r2;

        alpha2 = std::min(2 * alpha1, 0.5 * (alpha_max + alpha1));
        theta2 = Theta(spline, points, sw, alpha2, direction, fixed_knots);
        if (theta2 < 0)
            return false;

        q2 = mDelta;
        r2 = mPenalty;
    }

    // find Q coefficients
    double a0 = q0;
    double diff1 = alpha1 - alpha0;
    double diff2 = alpha2 - alpha0;
    double a2 = (q1 - q0) / diff1;
    a2 -= (q2 - q0) / diff2;
    a2 /= alpha1 - alpha2;
    double a1 = (q1 - a0) / diff1;
    a1 -= a2 * diff1;

    // find R coefficients
    double fraction = diff1 / diff2;
    double numerator = r1 - r0 - fraction * (r2 - r0);
    double temp = std::log((alpha_max - alpha1) / (alpha_max - alpha0));
    double denominator = temp - fraction * std::log((alpha_max - alpha2) / (alpha_max - alpha0));
    double b2 = numerator / denominator;
    double b1 = (r1 - r0 - b2 * temp) / diff1;

    // find coefficients of quadratic equation
    a = -2 * a2;
    b = -a * (alpha_max + alpha0) - mP * b1 - a1;
    double c = (mP * b1 + a1 + a * alpha0) * alpha_max - mP * b2;

    double root1 = -0.5 * (b + std::sqrt(b * b - 4 * a * c)) / a;
    double root2 = -b / a - root1;

    double alpha_res = 0;
    if (root1 > 0 && root1 < alpha_max)
        alpha_res = root1;
    else if (root2 > 0 && root2 < alpha_max)
        alpha_res = root2;

    double theta_res = Theta(spline, points, sw, alpha_res, direction, fixed_knots);

    if (theta_res < 0)
        return false;
    return true;
}

bool SplineController::InitiateGrid(Spline& spline, const SplinePoints& points, int nInternalKnots)
{
    int k = spline.GetDegree();
    int g = nInternalKnots;
    mKnots.setZero();
    mKnots[0] = points.front().x;
    mKnots[g + 1] = points.back().x;
    size_t n = points.size();

    int unique_size = 0;
    size_t index = 0;

    while (index < n && points[index].x < mKnots[g + 1])
    {
        if (index != 0 && !equal(points[index].x, points[index - 1].x))
            ++unique_size;
        ++index;
    }

    // not enough data points
    if (unique_size <= 0)
        return false;

    // number of knots should be less than n - k for n points with unique x
    if (unique_size < g + k + 1)
        return false;

    int points_per_knot = unique_size / (g + 1);
    int knot_index = 1;
    int i = 1;
    int counter = 0;

    while (knot_index < g + 1)
    {
        while (counter < knot_index * points_per_knot || equal(points[i].x, points[i - 1].x))
        {
            if (!equal(points[i].x, points[i - 1].x))
                ++counter;
            ++i;
        }
        mKnots[knot_index] = 0.5 * (points[i].x + points[i - 1].x);
        ++knot_index;
    }

    spline.SetKnots(mKnots);
    return true;
}

bool SplineController::ApproximateWithOptimalGrid(Spline& spline, const SplinePoints& points, double smooth_weight, double eps1, double eps2)
{
    int g = spline.GetInternalKnotsNum();
    int nk = mKnots.size();

    if (!Approximate(spline, points, smooth_weight))
        return false;

    Eigen::VectorXd direction(g);
    Eigen::VectorXd error_derivative(g);
    Delta(spline, points, smooth_weight);
    mP = eps1 * mDelta * (spline.GetRightBound() - spline.GetLeftBound()) / SQ(g + 1);
    mError = mDelta + mP * Penalty(spline);
    for (int i = 0; i < g; ++i)
    {
        error_derivative[i] = ErrorGradient(spline, points, smooth_weight, i + 1);
        direction[i] = -error_derivative[i];
    }

    double old_norm = direction.norm();
    double criteria1 = eps1 + eps2;
    double criteria2 = criteria1;
    int max_num_of_iter = 1000;

    int iteration = 0;
    double eps2_sq = SQ(eps2);

    Eigen::VectorXd fixed_knots(nk);
    Eigen::VectorXd knots(nk);
    while ((criteria1 >= eps1 || criteria2 >= eps2_sq) && iteration < max_num_of_iter)
    {
        spline.GetKnots(fixed_knots);
        double old_error = mError;

        // can't minimize further, knots are too close to each other
        if (!SpecDimensionalMinimization(spline, points, smooth_weight, direction, error_derivative, fixed_knots))
        {
            spline.SetKnots(fixed_knots);
            return Approximate(spline, points, smooth_weight);
        }

        for (int i = 0; i < g; ++i)
            error_derivative[i] = ErrorGradient(spline, points, smooth_weight, i + 1);

        double new_norm = error_derivative.norm();

        if (iteration % g == 0)
        {
            direction = -error_derivative;
        }
        else
        {
            direction *= (new_norm / old_norm);
            direction -= error_derivative;
        }

        spline.GetKnots(knots);

        double numerator = (knots - fixed_knots).segment(1, nk - 2).squaredNorm();
        double denominator = fixed_knots.segment(1, nk - 2).squaredNorm();

        criteria1 = std::abs(old_error - mError) / old_error;
        criteria2 = numerator / denominator;

        old_norm = new_norm;

        ++iteration;
    }

    return true;
}