#pragma once
#include "Base.h"
#include "IPlanarCurve.h"
#include <vector>

struct Alignment //평면선형
{
	enum class Type { Line, Arc, Clothoid, Parabola3D };
	Type type{ Type::Line };

	// Arc: k = 1 / R
	double Rx{ 0.0 };      // +오목, -볼록

	// Clothoid: k = s / A^2
	double A{ 0.0 };      // Clothoid Parameter
	double s0{ 0.0 };     // Clothoid Chain at Start

	// Parabola2D
	// - L: 곡선 길이(체인)
	// - G1, G2: 시작/끝 기울기 (dz/dx, 단위는 기울기)
	double L{ 0.0 };
	double G1{ 0.0 };
	double G2{ 0.0 };
};

struct LongitudinalSection //종단선형
{
	enum class Type { Line, Arc, Parabola2D };
	Type type{ Type::Line };

	// 기준점
	double x0{ 0.0 };   // 체인 기준 (Origin x)
	double z0{ 0.0 };   // 기준 z
	double sx0{ 0.0 };  // 기준점에서의 종단 기울기 dz/dx

	// Arc
	double Rv{ 0.0 }; // +: 오목, -: 볼록

	// Parabola2D
	double L{ 0.0 };  // 곡선 길이(체인 기준)
	double G1{ 0.0 }; // 곡선 시작 기울기
	double G2{ 0.0 }; // 곡선 끝 기울기

	double sxx() const // 상수 2차미분(곡률)
	{ 
		switch (type)
		{
		case Type::Line: return 0.0;
		case Type::Arc: return (Rv == 0.0 ? 0.0 : (1.0 / Rv));
		case Type::Parabola2D: return (L == 0.0 ? 0.0 : (G2 - G1) / L);
		}
		return 0.0;
	}
	double slope(double x) const // x에서의 기울기 sx(x)
	{
		const double dx = (x - x0);
		switch (type)
		{
		case Type::Line: return sx0;
		case Type::Arc:	return sx0 + sxx() * dx;
		case Type::Parabola2D: return G1 + sxx() * (dx);
		}
		return sx0;
	}
	double elevation(double x) const // x에서의 고도 z(x)
	{
		const double dx = (x - x0);
		switch (type)
		{
		case Type::Line: return z0 + sx0 * dx;
		case Type::Arc:
			// z = z0 + sx0*dx + 0.5*(1/Rv)*dx^2
			return z0 + sx0 * dx + 0.5 * sxx() * dx * dx;
		case Type::Parabola2D:
			// 포물선 수직곡선의 표준식
			// 시작점에서 G1, z'' = (G2-G1)/L 가 상수 → 동일한 2차식 형태
			// 참고: dx는 x - x0 (x0를 포물선 시작점으로 두는 것을 권장)
			return z0 + G1 * dx + 0.5 * sxx() * dx * dx;
		}
		return z0;
	}	
};

struct CrossSection // 횡단선형
{
	enum class Type { Line };
	Type type{ Type::Line };
};

class Surface
{
public:
	enum class Type 
	{
		Plane,                // only slopes, no curvature
		CurvatureNoTwist,     // has fxx, fyy but fxy = 0
		CurvatureWithTwist    // has fxx, fyy, fxy
	};

	Surface();
	virtual ~Surface();

	void Init();
	void RemoveAll();

	IPlanarCurve* pPlanarCrv;		// 평면선형
	LongitudinalSection* pLongi;	// 종단면
	CrossSection* pCross;			// 횡단면

	Point3d ptOrigin;		// (x0, y0, z0)
	double sx;				// Longitudinal slope (dz/dx)
	double sy;				// Transverse slope (dz/dy)
	double sxx;				// Longitudinal curvature (d2z/dx2)
	double syy;				// Transverse curvature (d2z/dy2)
	double sxy;				// Twist (d2z/dxdy)
	double tol{ 1e-12 };    // numerical tolerance for zero tests

public:
	void			ComputeSecondDerivsFromMeta();
	double			Elevation(double x, double y) const;
	Type			GetType() const;
	bool			IsPlane() const { return GetType() == Type::Plane; };
	bool			HasTwist() const { return std::abs(sxy) >= tol; };
	bool			HasAnyCurvature() const { return (std::abs(sxx) >= tol) || (std::abs(syy) >= tol); };

};