#include "Base.h"
#include "Arc.h"
#include <corecrt_math_defines.h>

Arc::Arc()
{
	Init();
}

Arc::~Arc()
{
	RemoveAll();
}

void Arc::Init()
{
	p0 = Point2d();
	th0 = 0.0;
	R = 0.0;
	_MaxSta = std::abs(R) * 2.0 * M_PI;
	_Sta = 0.0;
}

void Arc::RemoveAll()
{

}

static double wrapTo2Pi(double rad) 
{
	rad = std::fmod(rad, 2.0 * M_PI);
	if (rad < 0.0) rad += 2.0 * M_PI;

	return rad;
}

Point2d Arc::position(double s) const
{
	Point2d pt;

	double dth = s / R; //���ఢ
	Point2d C = GetCentroid();

	Vector2d v0 = Vector2d(- R * std::cos(th0), -R * std::sin(th0)); // C->p0
	Vector2d v = v0.rotate(dth); // ȸ��
	pt = Point2d(C.x + v.x, C.y + v.y); // �̵�

	return pt;
}

double Arc::theta(double s) const
{
	return th0 + s / R; //+����, -����
}

bool Arc::SetStaByPoint(const Point2d& q)
{
	// ���� �������� ����ǥ�� ����
	Point2d C = GetCentroid();
	Vector2d cq{ q.x - C.x, q.y - C.y };

	// ������ ���� üũ
	double r_q = std::hypot(cq.x, cq.y);
	double tol = 1e-6 * std::abs(R);
	if (std::abs(r_q - std::abs(R)) > tol)
		return false; // �� ���� ���� ����

	double ang_q = std::atan2(cq.y, cq.x);
	double ang_p0 = std::atan2(p0.y - C.y, p0.x - C.x);

	int dir = (R >= 0.0 ? +1 : -1);

	// ���ఢ (�׻� ���)
	double dth;
	if (dir > 0) 
	{
		dth = wrapTo2Pi(ang_q - ang_p0);
	}
	else 
	{
		dth = wrapTo2Pi(ang_p0 - ang_q);
	}

	double s = std::abs(R) * dth;
	if (s < -tol || s > _MaxSta + tol)
		return false;

	if (s < 0.0) s = 0.0;
	if (s > _MaxSta) s = _MaxSta;
	_Sta = s;

	return true;
}

Point2d Arc::GetCentroid() const
{
	Vector2d vP0ToCen{ -std::sin(th0), std::cos(th0) }; //������-> 90�� ȸ��

	return { p0.x + R * vP0ToCen.x, p0.y + R * vP0ToCen.y };
}