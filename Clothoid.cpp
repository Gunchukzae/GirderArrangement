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
	// ��ġ ����ȭ(Newton/��缽��)�� ||pos(s)-q|| �ּ� s �� [0,L]
	// ���⼭�� �������̽���: ���� ������ �ļ� �ܰ迡�� �߰�.
	return true;
}
