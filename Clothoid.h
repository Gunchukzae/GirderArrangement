#pragma once
#include "IPlanarCurve.h"
#include "Base.h"

// ---------------------------------------------
// Clothoid
// ---------------------------------------------
//   k(s) = s/A^2, θ(s) = th0 + s^2/(2A^2)
//   로컬좌표 Δp(s) = A*sqrt(2π) * ( [C(u)-C(0)], [S(u)-S(0)] ),  u = s/(sqrt(π)*A)
//   (C,S)는 Fresnel 적분

class Clothoid : public IPlanarCurve 
{
	Clothoid();
	virtual ~Clothoid();

	void Init();
	void RemoveAll();

public:
	Point2d position(double s) const override;
	double heading(double s) const override { return th0 + (s * s) / (2.0 * A * A); }
	double curvature(double s) const override {	return s / (A * A); }
	double length() const override { return L; }
	bool SetSByPoint(const Point2d& q) override;

public:
	static void fresnelCS(double u, double& C, double& S); // Fresnel 적분 근사(간단 placeholder; 실제로는 고정밀 근사 사용 권장)

public:
	Point2d p0;	// 시작점
	double th0;	// 시작 접선각(라디안)
	double A;	// 클로소이드 파라미터
	double L;	// 구간 길이
};