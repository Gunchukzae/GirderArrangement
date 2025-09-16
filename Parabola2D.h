#pragma once
#include "IPlanarCurve.h"
#include "Base.h"
#include <limits>

// ---------------------------------------------
// Parabola2D
// ---------------------------------------------
// 로컬좌표(ξ,η)에서  η = (1/(2p)) * ξ^2  (표준형; p>0는 초점거리)
// 시작점 P0, 시작 접선각 th0 기준으로 로컬-글로벌 변환
// 아크길이 s와 ξ 사이의 정확한 역변환은 닫힌형이 없어 수치해석 사용.

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
	// 로컬 좌표에서 파라메트릭 by ξ
	static inline Vector2d localByXi(double xi, double p) {
		double eta = (xi * xi) / (2.0 * p);
		return { xi, eta };
	}
	// 로컬 접선각 θ_local(xi) = atan(dη/dξ) = atan(xi/p)
	static inline double localHeadingByXi(double xi, double p) {
		return std::atan(xi / p);
	}
	// 로컬 곡률 k(xi) = |dθ/ds| = |(p)/( (p^2 + xi^2)^(3/2) )|
	static inline double localCurvatureByXi(double xi, double p) {
		return std::abs(p) / std::pow(p * p + xi * xi, 1.5);
	}
	// s(ξ) = ∫ sqrt(1 + (ξ/p)^2) dξ = (ξ/2)*sqrt(1+(ξ/p)^2) + (p/2) * asinh(ξ/p)
	static inline double arcLengthOfXi(double xi, double p) {
		double t = xi / p;
		return 0.5 * xi * std::sqrt(1.0 + t * t) + 0.5 * p * std::asinh(t);
	}
	// s -> ξ 수치역함수 (뉴턴/이분법)
	static double xiFromS(double s, double p);

private:
	Point2d m_P0;	 // 시작점
	double m_Theta0; // 시작 접선각(라디안)
	double m_Pos;	 // 초점거리 파라미터(>0)
	double m_Len;	 // 사용 구간의 체인 길이(0~L)
};
