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

	double dth = s / R; //진행각
	Point2d C = GetCentroid();

	Vector2d v0 = Vector2d(- R * std::cos(th0), -R * std::sin(th0)); // C->p0
	Vector2d v = v0.rotate(dth); // 회전
	pt = Point2d(C.x + v.x, C.y + v.y); // 이동

	return pt;
}

double Arc::theta(double s) const
{
	return th0 + s / R; //+오목, -볼록
}

bool Arc::SetStaByPoint(const Point2d& q)
{
	// 중점 기준으로 극좌표로 투영
	Point2d C = GetCentroid();
	Vector2d cq{ q.x - C.x, q.y - C.y };

	// 반지름 오차 체크
	double r_q = std::hypot(cq.x, cq.y);
	double tol = 1e-6 * std::abs(R);
	if (std::abs(r_q - std::abs(R)) > tol)
		return false; // 원 위에 있지 않음

	double ang_q = std::atan2(cq.y, cq.x);
	double ang_p0 = std::atan2(p0.y - C.y, p0.x - C.x);

	int dir = (R >= 0.0 ? +1 : -1);

	// 진행각 (항상 양수)
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
	Vector2d vP0ToCen{ -std::sin(th0), std::cos(th0) }; //접선각-> 90도 회전

	return { p0.x + R * vP0ToCen.x, p0.y + R * vP0ToCen.y };
}