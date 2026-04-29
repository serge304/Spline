/*
 * SplineController.h
 *
 *  Created on: Jul 25, 2022
 *      Author: sergey
 */

#ifndef SPLINECONTROLLER_H_
#define SPLINECONTROLLER_H_

#include "Spline.h"

class SplineController
{
public:
    SplineController();

private:
    double Penalty(const Spline& spline);

    double Delta(const Spline& spline, const SplinePoints& points, double smoothing_weight);

    double PenaltyDerivative(const Spline& spline, int knot_id);

    double Error(const Spline& spline, const SplinePoints& points, double sw);

    double ErrorGradient(const Spline& spline, const SplinePoints& points, double smoothing_weight, int knot_id);

    double Theta(Spline& spline, const SplinePoints& points, double sw, double alpha, const Eigen::VectorXd& direction, const Eigen::VectorXd& fixed_knots);

    bool SpecDimensionalMinimization(Spline& spline, const SplinePoints& points, double sw,
                                     const Eigen::VectorXd& direction, const Eigen::VectorXd& error_derivative,
                                     const Eigen::VectorXd& fixed_knots);

    bool InitiateGrid(Spline& spline, const SplinePoints& points, int nInternalKnots);

    void Reset();

    void SetBuffers(int nInternalKnots, int degree);

public:
    bool Init(Spline& spline, const SplinePoints& points);

    bool Init(Spline& spline, const SplinePoints& points, int nInternalKnots);

    bool Approximate(Spline& spline, const SplinePoints& points, double smoothing_weight);

    bool ApproximateWithOptimalGrid(Spline& spline, const SplinePoints& points, double smooth_weight, double eps1, double eps2);

private:
    double mDelta;
    double mP;
    double mError;
    double mPenalty;

    // буферы
    std::vector<double> mSplines;
    Eigen::VectorXd mKnots;
    Eigen::MatrixXd mA;
    Eigen::VectorXd mR, mC;
};

#endif /* SPLINECONTROLLER_H_ */