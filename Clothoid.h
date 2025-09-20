#pragma once
#include "IPlanarCurve.h"
#include "Base.h"

// ---------------------------------------------
// Clothoid (A)
// ---------------------------------------------

class Clothoid : public IPlanarCurve 
{
public:
	Clothoid();
	virtual ~Clothoid();

	void Init();
	void RemoveAll();

	Point2d position(double s) const override;
	double theta(double s) const override { return th0 + (s * s) / (2.0 * A * A); }
	double curvature(double s) const override {	return s / (A * A); }
	double length() const override { return _Sta; }
	bool SetStaByPoint(const Point2d& q) override;

private:
	void GetDistDeriveBetweenStaPt(double s, const Point2d& q, double& fp, double& fpp);
	Vector2d GetTanVecBySta(double s);
	Vector2d GetNormVecBySta(double s);
	double GetDist2BetweenStaPt(const double s, const Point2d& q);

public:
	Point2d p0;	// ������
	double th0;	// ���� ������(����)
	double A;	// Ŭ�μ��̵� �Ķ����

private:
	double _MaxSta; // Ŭ�μ��̵� Max Length
	double _Sta;	// ���� station (= p0->STA ����)
};