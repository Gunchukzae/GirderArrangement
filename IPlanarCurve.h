#pragma once
#include <cmath>

struct Point2d;
class IPlanarCurve 
{
public:
	virtual ~IPlanarCurve() = default;

	// station�� s�� ��
	virtual Point2d position(double s)   const = 0;	// ��ġ P(s)
	virtual double heading(double s) const = 0;		// ������ ��(s)
	virtual double curvature(double s) const = 0;	// ��� k(s) = d��/ds
	virtual double length() const = 0;				// ������� l(s)
	virtual bool SetSByPoint(const Point2d& q) = 0; // q�� station ����

};