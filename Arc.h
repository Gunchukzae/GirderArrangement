#pragma once
#include "IPlanarCurve.h"
#include "Base.h"

// ---------------------------------------------
// Arc (R, +: CCW, -: CW)
// ---------------------------------------------

class Arc : public IPlanarCurve
{
public:
	Arc();
	virtual ~Arc();

	void Init();
	void RemoveAll();

	Point2d position(double s) const override;
	double theta(double s) const override;
	double curvature(double /*s*/) const override { return 1.0 / R;	}
	double length() const override { return _Sta; }
	bool SetSByPoint(const Point2d& q) override;

	Point2d GetCentroid() const;

public:
	Point2d p0;	// 시작점
	double th0;	// 시작 접선각(라디안)
	double R;	// R (+:CCW, -:CW)

private:
	double  _MaxSta; // 호 Max Length
	double  _Sta;	 // 현재 station (= p0->STA 길이)
};