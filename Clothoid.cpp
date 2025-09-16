#include "Clothoid.h"
#include "Base.h"
#include <corecrt_math_defines.h>
#include <limits>

Clothoid::Clothoid()
{
	Init();
}

Clothoid::~Clothoid()
{
	RemoveAll();
}

void Clothoid::Init()
{
	p0 = Point2d();
	th0 = 0.0;
	A = 0.0;
	L = 0.0;
}

void Clothoid::RemoveAll()
{

}

Point2d Clothoid::position(double s) const
{
	double u = s / (std::sqrt(M_PI) * A);
	double C, S, C0, S0;
	fresnelCS(u, C, S);
	fresnelCS(0.0, C0, S0); // = (0,0)
	double scale = A * std::sqrt(2.0 * M_PI);
	Vector2d dL{ scale * (C - C0), scale * (S - S0) };

	Vector2d dG = dL.rotate(th0);
	return { p0.x + dG.x, p0.y + dG.y };
}

bool Clothoid::SetSByPoint(const Point2d& q)
{
	// 수치 최적화(Newton/골든섹션)로 ||pos(s)-q|| 최소 s ∈ [0,L]
	// 여기서는 인터페이스만: 실제 구현은 후속 단계에서 추가.
	return true;
}
