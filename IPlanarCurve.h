#pragma once
#include <cmath>

struct Point2d;
class IPlanarCurve 
{
public:
	virtual ~IPlanarCurve() = default;

	// station이 s일 때
	virtual Point2d position(double s)   const = 0;	// 위치 P(s)
	virtual double heading(double s) const = 0;		// 접선각 θ(s)
	virtual double curvature(double s) const = 0;	// 곡률 k(s) = dθ/ds
	virtual double length() const = 0;				// 진행길이 l(s)
	virtual bool SetSByPoint(const Point2d& q) = 0; // q점 station 설정

};