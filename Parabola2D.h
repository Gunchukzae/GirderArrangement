#pragma once
#include "IPlanarCurve.h"
#include "Base.h"
#include <limits>

// ---------------------------------------------
// Parabola2D
// ---------------------------------------------
// ������ǥ(��,��)����  �� = (1/(2p)) * ��^2  (ǥ����; p>0�� �����Ÿ�)
// ������ P0, ���� ������ th0 �������� ����-�۷ι� ��ȯ
// ��ũ���� s�� �� ������ ��Ȯ�� ����ȯ�� �������� ���� ��ġ�ؼ� ���.

class Parabola2D : public IPlanarCurve
{
	Parabola2D();
	virtual ~Parabola2D();

	void Init();
	void RemoveAll();

public:
	Point2d position(double s) const override;
	double heading(double s) const override;
	double curvature(double s) const override;
	double length() const override { return m_Len; }
	bool SetSByPoint(const Point2d& q) override;

private:
	// ���� ��ǥ���� �Ķ��Ʈ�� by ��
	static inline Vector2d localByXi(double xi, double p) {
		double eta = (xi * xi) / (2.0 * p);
		return { xi, eta };
	}
	// ���� ������ ��_local(xi) = atan(d��/d��) = atan(xi/p)
	static inline double localHeadingByXi(double xi, double p) {
		return std::atan(xi / p);
	}
	// ���� ��� k(xi) = |d��/ds| = |(p)/( (p^2 + xi^2)^(3/2) )|
	static inline double localCurvatureByXi(double xi, double p) {
		return std::abs(p) / std::pow(p * p + xi * xi, 1.5);
	}
	// s(��) = �� sqrt(1 + (��/p)^2) d�� = (��/2)*sqrt(1+(��/p)^2) + (p/2) * asinh(��/p)
	static inline double arcLengthOfXi(double xi, double p) {
		double t = xi / p;
		return 0.5 * xi * std::sqrt(1.0 + t * t) + 0.5 * p * std::asinh(t);
	}
	// s -> �� ��ġ���Լ� (����/�̺й�)
	static double xiFromS(double s, double p);

private:
	Point2d m_P0;	 // ������
	double m_Theta0; // ���� ������(����)
	double m_Pos;	 // �����Ÿ� �Ķ����(>0)
	double m_Len;	 // ��� ������ ü�� ����(0~L)
};
