#pragma once
#include "IPlanarCurve.h"
#include "Base.h"
#include <limits>

// ---------------------------------------------
// Parabola2D (+P(∪)/-P(∩)): η = ξ^2/2P (ξ=xi, η=eta)
// ---------------------------------------------

class Parabola2D : public IPlanarCurve
{
public:
	Parabola2D();
	virtual ~Parabola2D();

	void Init();
	void RemoveAll();

	Point2d position(double s) const override;
	double theta(double s) const override;
	double curvature(double s) const override;
	double length() const override { return _Sta; }
	bool SetStaByPoint(const Point2d& q) override;

private:
	double GetDist2BetweenStaPt(double s, const Point2d& q) const;
	Vector2d GetTanVecBySta(double s);
	Vector2d GetNormVecBySta(double s);
	void GetDistDeriveBetweenStaPt(double s, const Point2d& q, double& fp, double& fpp);
	double EstimateArcLengthUpperBound(const Point2d& q);

public:
	Point2d p0;		// 시작점
	double th0;		// 시작 접선각(라디안)
	double P;		// 초점거리 파라미터(>0)

private:
	double _Sta;	// 현재 station (= p0->STA 길이)
};
