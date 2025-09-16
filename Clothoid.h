#pragma once
#include "IPlanarCurve.h"
#include "Base.h"

// ---------------------------------------------
// Clothoid
// ---------------------------------------------
//   k(s) = s/A^2, ��(s) = th0 + s^2/(2A^2)
//   ������ǥ ��p(s) = A*sqrt(2��) * ( [C(u)-C(0)], [S(u)-S(0)] ),  u = s/(sqrt(��)*A)
//   (C,S)�� Fresnel ����

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
	static void fresnelCS(double u, double& C, double& S); // Fresnel ���� �ٻ�(���� placeholder; �����δ� ������ �ٻ� ��� ����)

public:
	Point2d p0;	// ������
	double th0;	// ���� ������(����)
	double A;	// Ŭ�μ��̵� �Ķ����
	double L;	// ���� ����
};