#include "CurveFunctions.h"
using std::vector;
#include <iostream>  
#include <sstream>
#include <Eigen/LU>
using Eigen::MatrixXd;
using Eigen::Matrix4d;
using Eigen::Vector4d;

#include <cassert>
#include <cmath>
#include <math.h> 

#include "jsassert.h"
#undef assert
#define assert(cond) jsAssert(cond)

// Call these to raise a dialog box or log to the javascript console for debugging.
// NOTE: You can pass either a const char* or an std::string.
// void jsAlert( const std::string& msg );
// void jsLog( const std::string& msg );
// void jsWarn( const std::string& msg );
// void jsError( const std::string& msg );

namespace Curve
{

// Bezier helper functions.
namespace
{
// You may add helper functions here.



}


// Evaluate a cubic Bezier curve at location 't'.
Point EvaluateCubicBezierCurve( const Point& p0, const Point& p1, const Point& p2, const Point& p3, const real_t t, EvaluateCubicBezierCurveApproach approach )
{
    if( BernsteinApproach == approach ) return EvaluateCubicBezierCurveBernstein( p0, p1, p2, p3, t );
    else if( MatrixApproach == approach ) return EvaluateCubicBezierCurveMatrix( p0, p1, p2, p3, t );
    else if( CasteljauApproach == approach ) return EvaluateCubicBezierCurveCasteljau( p0, p1, p2, p3, t );
    else {
        assert( !"Unknown EvaluateCubicBezierCurveApproach" );
        return Point(-31337,-31337);
    }
}
Point EvaluateCubicBezierCurveBernstein( const Point& p0, const Point& p1, const Point& p2, const Point& p3, const real_t t )
{
    
	Point p = pow(1.0 - t, 3.0) * (p0) + 3.0 * t*pow(1.0 - t, 2.0)*p1 + 3.0* (pow(t, 2.0))*(1.0 - t)*p2 + (pow(t, 3.0))*p3;

	return p;
}
Point EvaluateCubicBezierCurveMatrix( const Point& p0, const Point& p1, const Point& p2, const Point& p3, const real_t t )
{
  
	MatrixXd A(4, 4);
	MatrixXd p(1, 4);

	MatrixXd prod(1,4);


	// basis matrix
	A(0, 0) = -1.0; A(0, 1) = 3.0; A(0, 2) = -3.0; A(0, 3) = 1.0;
	A(1, 0) = 3.0; A(1, 1) = -6.0; A(1, 2) = 3.0; A(1, 3) = 0.0;
	A(2, 0) = -3.0; A(2, 1) = 3.0; A(2, 2) = 0.0; A(2, 3) = 0.0;
	A(3, 0) = 1.0; A(3, 1) = 0.0; A(3, 2) = 0.0; A(3, 3) = 0.0;

	// t matrix
	p(0, 0) = pow(1.0-t, 3.0); p(0, 1) = pow(1.0-t, 2.0); p(0, 2) = 1.0-t; p(0, 3) = 1;
	prod = p * A;



	Point r;
	r(0) = prod(0, 0)*p3(0) + prod(0,1)*p2(0) + prod(0,2)*p1(0) + prod(0,3)*p0(0);
	r(1) = prod(0, 0)*p3(1) + prod(0, 1)*p2(1) + prod(0, 2)*p1(1) + prod(0, 3)*p0(1);
	

    return r;
}
Point EvaluateCubicBezierCurveCasteljau( const Point& p0, const Point& p1, const Point& p2, const Point& p3, const real_t t )
{
    
	Point r;
	// as directed in slides
	r = (1 - t)*((1 - t)*((1 - t)*p0 + t * p1) + t * ((1 - t)*p1 + t * p2)) + t * ((1 - t)*((1 - t)*p1 + t * p2) + t * ((1 - t)*p2 + t * p3));
    return r;
}

// Evaluate a cubic Bezier spline with control points 'controlPoints' arranged
//     on_curve ( off_curve off_curve on_curve )+
// at positive integer 'samplesPerCurve' locations along each curve.
// Returns the sampled points.
std::vector< Point > EvaluateCubicBezierSpline( const std::vector< Point >& controlPoints, const int samplesPerCurve, EvaluateCubicBezierCurveApproach approach )
{
    assert( controlPoints.size() >= 4 );
    assert( samplesPerCurve > 0 );

	//const std::vector< Point >& C = controlPoints;
	
	std::vector< Point > result;
	
	int n = controlPoints.size();
	int j = 0;
	//std::stringstream ss(std::stringstream::in | std::stringstream::out);
	//const std::string& msg = "this runs";
	
		for (j = 0; j  < n-1; j+=3) {
			real_t t = 0.0;
			for (int i = 0; i != samplesPerCurve; i++) {
				
			


				result.push_back(EvaluateCubicBezierCurve(controlPoints[j], controlPoints[j + 1], controlPoints[j + 2], controlPoints[j + 3], t, approach));
				t = t + (1.0 / samplesPerCurve);
			
			//ss << j;
			//const std::string& test = ss.str();
			//jsAlert(test);

		}

	}
    // ADD YOUR CODE HERE
   
    return result;
}

/// ======================================================================================

// Evaluate a cubic Hermite curve at location 't'.
Point EvaluateCubicHermiteCurve( const Point& p0, const Point& dp0, const Point& p1, const Point& dp1, const real_t t )
{
    // ADD YOUR CODE HERE
	Point r;
	

	MatrixXd A(4, 4);
	MatrixXd p(4, 1);
	MatrixXd h(2, 4);

	MatrixXd prod(4, 1);

	

	A(0, 0) = -2.0; A(0, 1) = 3.0; A(0, 2) = 0.0; A(0, 3) = 0.0;
	A(1, 0) = 2.0; A(1, 1) = -3.0; A(1, 2) = 0.0; A(1, 3) = 1.0;
	A(2, 0) = 1.0; A(2, 1) = -1.0; A(2, 2) = 0.0; A(2, 3) = 0.0;
	A(3, 0) = 1.0; A(3, 1) = -2.0; A(3, 2) = 1.0; A(3, 3) = 0.0;

	p(0, 0) = pow(t, 3.0); p(1, 0) = pow(t, 2.0); p(2,0) = t; p(3,0) = 1;
	prod = A*p;

	h(0, 0) = p1(0,0); h(0, 1) = p0(0,0); h(0, 2) = dp1(0,0); h(0, 3) = dp0(0,0);
	h(1, 0) = p1(1,0); h(1, 1) = p0(1,0); h(1, 2) = dp1(1,0); h(1, 3) = dp0(1,0);

	

	r(0) = h(0, 0)*prod(0, 0) + h(0, 1)*prod(1, 0) + h(0, 2)*prod(2, 0) + h(0, 3)*prod(3, 0);
	r(1) = h(1, 0)*prod(0, 0) + h(1, 1)*prod(1, 0) + h(1, 2)*prod(2, 0) + h(1, 3)*prod(3, 0);
	

	//-------------------------------------------------------
	/*
	 another approach without matrix multiplication
	float h1= 2*pow(t,3.0) - 3*pow(t,2.0) + 1;
	float h2= -2*pow(t,3.0) + 3*pow(t,2.0);
	float h3 = pow(t,3.0) - 2 * pow(t,2.0) + t;
	float h4 = pow(t,3.0) - pow(t,2.0);

	r = h1 * p0 + h2 * p1 + h3 * dp0 + h4 * dp1;
	*/

	//-------------------------------------------------------

    return r;
}

// Evaluate a cubic Hermite spline with control points 'controlPoints' arranged:
//     p0 derivative_at_p0 ( p1 derivative_at_p1 )+
// at positive integer 'samplesPerCurve' locations along each curve.
// Upon return, 'curvePointsOut' is cleared and replaced with the sampled points.
std::vector< Point > EvaluateCubicHermiteSpline( const std::vector< Point >& controlPoints, const int samplesPerCurve )
{
    assert( controlPoints.size() >= 4 );
    assert( samplesPerCurve > 0 );
    
    // ADD YOUR CODE HERE
    std::vector< Point > result;

	int j = 0;

	//const std::string& msg = "this works";
	
	
		while (j + 3 < controlPoints.size()) {


			double t = 0.0;

				for (int i = 0; i <= samplesPerCurve; i++) {
					
					result.push_back(EvaluateCubicHermiteCurve(controlPoints[j], controlPoints[j + 1], controlPoints[j + 2], controlPoints[j + 3], t));
					t = t + (1.0 / samplesPerCurve);
				}
					
				
			j=j+2;
			}
			//jsAlert(msg);

	
    return result;
}

// Given a cubic Hermite spline with control points 'controlPoints' arranged:
//     p0 derivative_at_p0 ( p1 derivative_at_p1 )+
// replaces the derivative entries with values that result in a C2 continuous
// Hermite spline.
// NOTE: 'controlPoints' is an input and output parameter. The derivative entries are replaced.

// Hint: To implement C2 continuous Hermite splines, you need to solve a system of equations.
//       This is best done with a matrix.
//       Then you solve a linear system Ac = p, where c are the unknown derivatives.
//       I have included the Eigen linear algebra package that can solve linear systems.
//       Below is an example showing how to use the linear system solver.
//
//  MatrixXd A(3,3);
//  MatrixXd c(3,1);
//  MatrixXd p(3,1);
//  A(0,0) = 1.0; A(0,1) = 0.0; A(0,2) = 0.0;
//  A(1,0) = 0.0; A(1,1) = 1.0; A(1,2) = 0.0;
//  A(2,0) = 0.0; A(2,1) = 0.0; A(2,2) = 1.0;
//  p(0,0) = 1.0; p(1,0) = 2.0; p(2,0) = 3.0;
// c = A.fullPivLu().solve(p);
//
//  The result will be stored in c as follows: c(0,0) = 1.0; c(1,0) = 2.0; c(3,0) = 3.0, which satisfies Ac = p.
void CalculateHermiteSplineDerivativesForC2Continuity( std::vector< Point >& controlPoints )
{
    // Do nothing if there aren't enough control points.
    if( controlPoints.size() < 4 ) return;
    
    assert( controlPoints.size() >= 4 );
    assert( controlPoints.size() % 2 == 0 );

	//-------------------------------------------------START------------------------------------------------------------
	//remove below /* to make this function work and the */ before the END line

	/*

	//the idea is to generate a basis matrix with the coeeficients generated by the equations and the p matrix with points according to the equaltion given
	int N = controlPoints.size();
	
	MatrixXd A(N, N);
	MatrixXd p(N, 2);
	MatrixXd c(N, 2);
	//the A matrix (basis matrix) has the values from i=1 to i= N-2 like 1 4 1 according to what 'd'  it will be multiplied to, like the equations are d0 + 4d1 +d2 = 3*(p2-p0)
	for (int i = 0; i <= N - 1; i++) {
		for (int j = 0; j <= N-1; j++) {
			A(i, j) = 0.0;
		}
	}



	//the first and the last row of the matrix A are hardcoded with values 1 2 .... as per the second derivative equations of first piece and last piece equating to 0
	A(0, 0) = 2.0; A(0, 1) = 1.0;
	A(N - 2, N - 3) = 2.0; A(N - 2, N - 2) = 1.0;
	for (int i = 1; i <= N - 2; i++) {
	
	
			A(i, i-1) = 1.0;
			A(i, i ) = 4.0;
			A(i, i + 1) = 1.0;

	}

	the p matrix is also hardcoded for first and last row and is filled in by the loop for in between first and last row
	 p.row(0) = 3 * (controlPoints[2] - controlPoints[0]);

	int j = 2;

	

	for (int i = 0; i <= N - 2; i+=2) {
		p.row(i) = 3 * (controlPoints[i+2] - controlPoints[i]);

	j = j + 2;
	}

	p.row(N-1) = 3 * (controlPoints[N-2] - controlPoints[N-4]);

	// solve for c and assign all the derivatives to controlPoints input vector for i= 1,3,5 ... N-1
	c = A.fullPivLu().solve(p);

	int k = 0;
	for (int x = 1; x <= N - 1; x += 2) {
		controlPoints[x] = c.row(k);
		
		k++;
	}
	
	*/


	//-------------------------------------------------END------------------------------------------------------------

}

/// ======================================================================================

// Evaluate a Catmull-Rom Spline with control points 'controlPoints' arranged:
//     p0 p1 p2 ( p3 )+
// at positive integer 'samplesPerCurve' locations along each curve.
// Upon return, 'curvePointsOut' is cleared and replaced with the sampled points.
std::vector< Point > EvaluateCatmullRomSpline(const std::vector< Point >& controlPoints, const int samplesPerCurve, const real_t alpha)
{
	assert(controlPoints.size() >= 4);
	assert(samplesPerCurve > 0);

	// ADD YOUR CODE HERE
	std::vector< Point > result;


	int j = 0;

	//const std::string& msg = "this works";

	//to make sure curve starts from first point, not working
	double t1 = 0.0;

	for (int i = 0; i <= samplesPerCurve; i++) {
		
		result.push_back(EvaluateCatmullRomCurve(2 * controlPoints[0] - controlPoints[1], controlPoints[0], controlPoints[1], controlPoints[2], t1, alpha));
		t1 = t1 + (1.0 / samplesPerCurve);
	}
	
	while (j + 3 < controlPoints.size()) {

		double t = 0.0;

			for (int i = 0; i <= samplesPerCurve; i++) {
				t = t + (1.0 / samplesPerCurve);
				result.push_back(EvaluateCatmullRomCurve(controlPoints[j], controlPoints[j + 1], controlPoints[j + 2], controlPoints[j + 3],t, alpha));

			}
			
		j++;
	}

	int k = controlPoints.size();
	double t2 = 0.0;
	for (int i = 0; i <= samplesPerCurve; i++) {

		result.push_back(EvaluateCatmullRomCurve(controlPoints[k - 3], controlPoints[k - 2], controlPoints[k - 1], 2 * controlPoints[k - 1] - controlPoints[k - 2], t2, alpha));
		t2 = t2 + (1.0 / samplesPerCurve);

	}

    return result;
}
// Evaluate a cubic Catmull-Rom Spline curve at location 't'.
Point EvaluateCatmullRomCurve( const Point& p0, const Point& p1, const Point& p2, const Point& p3, const real_t t, const real_t alpha )
{

	Point r;
	//-------------------------------------------------------------------------
	
	MatrixXd A(4, 4);
	MatrixXd p(1, 4);
	MatrixXd h(4, 2);
	MatrixXd prod(1, 4);



	A(0, 0) = 0.0; A(0, 1) = 2.0; A(0, 2) = 0.0; A(0, 3) = 0.0;
	A(1, 0) = -1.0; A(1, 1) = 0.0; A(1, 2) = 1.0; A(1, 3) = 0.0;
	A(2, 0) = 2.0; A(2, 1) = -5.0; A(2, 2) = 4.0; A(2, 3) = -1.0;
	A(3, 0) = -1.0; A(3, 1) = 3.0; A(3, 2) = -3.0; A(3, 3) = 1.0;


	p(0, 0) = 1.0; p(0, 1) = t; p(0, 2) = pow(t, 2.0); p(0, 3) = pow(t, 3.0);

	prod = p * A;
	r = alpha * (prod(0, 0)*p0 + prod(0, 1)*p1 + prod(0, 2)*p2 + prod(0, 3)*p3);


	//--------------------------------------------direct implementation-------------------------
	//r = alpha *((2 * p1) +(-p0 + p2) * t + (2 * p0 - 5 * p1 + 4 * p2 - p3) * pow(t,2.0) + (-p0 + 3 * p1 - 3 * p2 + p3) * pow(t, 3.0));

	//-----------------------------------implementation of formula given in notes, slow in running---------------------------------

	// curve not proper
	/*double d1 = pow(pow(p1(0, 0) - p0(0, 0), 2.0) + pow(p1(1, 0) - p0(1, 0), 2.0), 0.5);
	double d2 = pow(pow(p2(0, 0) - p1(0, 0), 2.0) + pow(p2(1, 0) - p1(1, 0), 2.0), 0.5);
	double d3 = pow(pow(p3(0, 0) - p2(0, 0), 2.0) + pow(p3(1, 0) - p2(1, 0), 2.0), 0.5);

	Point B0 = p1;
	Point B1 = ((pow(d1, 2.0*alpha)*p2) - (pow(d2, 2.0*alpha)*p0) + (2*(pow(d1, 2.0*alpha) + 3 * pow(d1, alpha)*pow(d2, alpha) + pow(d2, 2.0*alpha))*p1)) / (3 * pow(d1, alpha)*(pow(d1, alpha) + pow(d2, alpha)));
	Point B2 = ((pow(d3, 2.0*alpha)*p1) - (pow(d2, 2.0*alpha)*p3) + (2*(pow(d3, 2.0*alpha) + 3 * pow(d3, alpha)*pow(d2, alpha) + pow(d2, 2.0*alpha))*p2)) / (3 * pow(d3, alpha)*(pow(d3, alpha) + pow(d2, alpha)));
	Point B3 = p2;

	r = B0*pow(1.0 - t, 3.0) + B1* t*pow(1.0 - t, 2.0) + B2* (pow(t, 2.0))*(1.0 - t) + B3*(pow(t, 3.0));
	*/
	//----------------------------------------------------------------------------------------------------------------------

    return r;
}

/// ======================================================================================

// B-Spline helper functions
namespace
{
// You may add helper functions here.
}

// Evaluate a cubic B-Spline with control points 'controlPoints' arranged:
//     p0 p1 p2 ( p3 )+
// at positive integer 'samplesPerCurve' locations along each curve.
// Upon return, 'curvePointsOut' is cleared and replaced with the sampled points.
std::vector< Point > EvaluateCubicBSpline( const std::vector< Point >& controlPoints, const int samplesPerCurve )
{
    assert( controlPoints.size() >= 4 );
    assert( samplesPerCurve > 0 );
    
    // ADD YOUR CODE HERE
    std::vector< Point > result;
    return result;
}
// Evaluate a cubic B-Spline curve at location 't'.
Point EvaluateCubicBSplineCurve( const Point& p0, const Point& p1, const Point& p2, const Point& p3, const real_t t )
{
    // ADD YOUR CODE HERE
    return Point(0,0);
}

// Compute cubic BSpline control points that interpolate the given points.
std::vector< Point > ComputeBSplineFromInterpolatingPoints( const std::vector< Point >& interpPoints )
{
    assert( interpPoints.size() >= 2 );
    
    // ADD YOUR CODE HERE
    std::vector< Point > controlPoints;
    return controlPoints;
}
// Given a sequence of cubic BSpline control points, returns the interpolating points
// that could have been used to create them via ComputeBSplineFromInterpolatingPoints().
std::vector< Point > ComputeInterpolatingPointsFromBSpline( const std::vector< Point >& controlPoints )
{
    assert( controlPoints.size() >= 4 );
    const std::vector< Point >& C = controlPoints;
    
    std::vector< Point > result;
    // The interpolated points are at the start and end of each cubic B-Spline
    // (they are continuous).
    // So let's just sample the t=1 point on every cubic BSpline,
    // as well as the t=0 point of the first one.
    result.push_back( EvaluateCubicBSplineCurve( C[0], C[1], C[2], C[3], 0. ) );
    for( int i = 0; i+3 < C.size(); ++i )
    {
        result.push_back( EvaluateCubicBSplineCurve( C[i], C[i+1], C[i+2], C[i+3], 1. ) );
    }
    
    return result;
}

}
