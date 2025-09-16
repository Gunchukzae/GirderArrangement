#include "Parabola2D.h"

Parabola2D::Parabola2D()
{
	Init();
}

Parabola2D::~Parabola2D()
{
	RemoveAll();
}

void Parabola2D::Init()
{
	m_P0 = Point2d();
	m_Theta0 = 0.0;
	m_Pos = 1.0;
	m_Len = 0.0;
}

void Parabola2D::RemoveAll()
{

}

Point2d Parabola2D::position(double s) const
{
	double xi = xiFromS(s, m_Pos);
	Vector2d Lc = localByXi(xi, m_Pos);

	Vector2d G = Lc.rotate(m_Theta0);
	return { m_P0.x + G.x, m_P0.y + G.y };
}

double Parabola2D::heading(double s) const
{
	double xi = xiFromS(s, m_Pos);
	return m_Theta0 + localHeadingByXi(xi, m_Pos);
}

double Parabola2D::curvature(double s) const
{
	double xi = xiFromS(s, m_Pos);
	return localCurvatureByXi(xi, m_Pos);
}

bool Parabola2D::SetSByPoint(const Point2d& q)
{
	// 수치 최적화로 구현(후속 단계). 여기서는 인터페이스만.
	return true;
}
