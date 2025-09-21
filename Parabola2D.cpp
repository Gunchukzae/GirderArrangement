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
	p0 = Point2d();
	th0 = 0.0;
	P = 1.0;
	_Sta = 0.0;
}

void Parabola2D::RemoveAll()
{

}

static Vector2d localByXi(double xi, double p) 
{
	Vector2d v;
	double eta = (xi * xi) / (2.0 * p);

	v = Vector2d(xi, eta);

	return v;
}

static double localThetaByXi(double xi, double p) 
{
	return std::atan(xi / p); // 기울기(xi) = xi/p
}

static double localCurvatureByXi(double xi, double p) 
{
	// k(xi) = |dθ/ds|
	// dθ/dξ = p / (p^2 + ξ^2)
	// ds/dξ = sqrt(1 + (ξ/p)^2)
	return std::abs(p) / std::pow(p * p + xi * xi, 1.5);
}

static double arcLengthOfXi(double xi, double p) 
{
	// s(ξ) = ∫ sqrt(1 + (ξ/p)^2) dξ
	double t = xi / p;

	return 0.5 * xi * std::sqrt(1.0 + t * t) + 0.5 * p * std::asinh(t);
}

static double xiFromS(double s, double p)
{
	if (s <= 0) return 0.0;

	const double absP = std::abs(p);
	auto S = [&](double xi) { return arcLengthOfXi(xi, p); };
	auto dS = [&](double xi) { return std::sqrt(1.0 + (xi / p) * (xi / p)); };

	// xi 초기화 (P값 기준 -> 작은 값 or 큰 값)
	double xi = (s < absP) ? s : std::sqrt(2.0 * absP * s);
	// [min,max] 설정
	double min = 0.0, max = std::max(xi * 2.0, std::sqrt(2.0 * absP * s) + 2.0 * absP);
	//(hi > s)
	while (S(max) < s) max *= 2.0;

	for (int it = 0; it < 20; ++it) 
	{
		// Newton-Raphson Method
		double f = S(xi) - s;
		double fp = dS(xi);

		double xi_new = xi - f / fp;

		// 브래킷 유지 & 갱신
		if (f > 0) max = xi;
		else       min = xi;

		// 브래킷 밖이면 이분법
		if (xi_new < min || xi_new > max)
			xi_new = 0.5 * (min + max);

		if (std::abs(xi_new - xi) < 1e-12) 
			return xi_new;

		xi = xi_new;
	}

	return xi; // return 안 했을 경우
}

Point2d Parabola2D::position(double s) const
{
	Point2d pt;

	double xi = xiFromS(s, P);
	Vector2d Lc = localByXi(xi, P);
	Vector2d G = Lc.rotate(th0); // 회전
	pt = Point2d(p0.x + G.x, p0.y + G.y); // 이동

	return pt;
}

double Parabola2D::theta(double s) const
{
	double xi = xiFromS(s, P);
	
	return th0 + localThetaByXi(xi, P);
}

double Parabola2D::curvature(double s) const
{
	double xi = xiFromS(s, P);
	
	return localCurvatureByXi(xi, P);
}

static double Clamp(double x, double a, double b)
{
	return x < a ? a : (x > b ? b : x);
}

Vector2d Parabola2D::GetTanVecBySta(double s)
{
	double th = theta(s);

	return Vector2d(std::cos(th), std::sin(th));
}

Vector2d Parabola2D::GetNormVecBySta(double s)
{
	double th = theta(s);

	return Vector2d(-std::sin(th), std::cos(th));
}

double Parabola2D::GetDist2BetweenStaPt(double s, const Point2d& q) const
{
	Point2d p = position(s);
	double dx = p.x - q.x;
	double dy = p.y - q.y;

	return dx * dx + dy * dy;
}

void Parabola2D::GetDistDeriveBetweenStaPt(double s, const Point2d& q, double& fp, double& fpp)
{
	Point2d p = position(s);
	double th = theta(s);
	double k = curvature(s);

	Vector2d r = Vector2d(p.x - q.x, p.y - q.y); // r(s) = p(s) − q
	Vector2d T = GetTanVecBySta(s);
	Vector2d N = GetNormVecBySta(s);

	// f(s) = r^2
	// f'(s) = 2 r·T
	fp = 2.0 * (r.x * T.x + r.y * T.y);

	// f''(s) = 2( |T|^2 + r·(k N) ) = 2(1 + k (r·N))
	fpp = 2.0 * (1.0 + k * (r.x * N.x + r.y * N.y));
}

double Parabola2D::EstimateArcLengthUpperBound(const Point2d& q)
{
	// [0,L] -> 최대 L값 추정 (포물선용 휴리스틱)
	const double _p = std::max(1e-12, P);
	const double dx = q.x - p0.x;
	const double dy = q.y - p0.y;
	const double D = std::sqrt(dx * dx + dy * dy);

	double xi_max = std::sqrt(std::max(0.0, 2.0 * _p * std::max(1.0, D)));
	xi_max *= 2.0; // 보정
	double L = arcLengthOfXi(xi_max, _p);

	// 안정성용 하한값 설정
	if (!std::isfinite(L) || L <= 0.0) L = 1.0;

	return L;
}

bool Parabola2D::SetStaByPoint(const Point2d& q)
{
	// 탐색 구간 길이 L 설정 (무한곡선이라 휴리스틱)
	const double L = EstimateArcLengthUpperBound(q);
	if (L <= 0.0 || !std::isfinite(L)) { _Sta = 0.0; return false; }

	// s 초기값 설정 (= 시작점 접선 투영길이)
	Vector2d T0 = Vector2d(std::cos(th0), std::sin(th0));
	Vector2d r0 = Vector2d(q.x - p0.x, q.y - p0.y);
	double s = Clamp(r0.x * T0.x + r0.y * T0.y, 0.0, L);

	// 1.Newton-Raphson Method
	const int    kMaxNewton = 12;
	const double dTolS = 1e-10 * std::max(1.0, L);
	const double dTolGrad = 1e-10;

	bool bSuccess = false;
	for (int it = 0; it < kMaxNewton; ++it)
	{
		double fp, fpp;
		GetDistDeriveBetweenStaPt(s, q, fp, fpp);

		if (std::fabs(fp) < dTolGrad)
		{
			bSuccess = true;
			break;
		}
		
		if (!(fpp > 0.0) || !std::isfinite(fpp)) break; // Error

		double step = fp / fpp;
		double s_new = s - step;

		// [0,L] 벗어났을 경우
		if (s_new < 0.0 || s_new > L)
		{
			double fCur = GetDist2BetweenStaPt(s, q);
			double f0 = GetDist2BetweenStaPt(0.0, q);
			double fL = GetDist2BetweenStaPt(L, q);

			// [0,L]에서 0,L이 최소 후보일 경우
			if (f0 < fCur || fL < fCur) 
			{
				s = (f0 < fL) ? 0.0 : L;
				bSuccess = true;
				break;
			}
			// 다시 시도
			s_new = Clamp(s_new, 0.0, L);
		}

		if (std::fabs(s_new - s) < dTolS)
		{
			s = s_new;
			bSuccess = true;
			break;
		}

		s = s_new;
	}

	// 2.Golden-section search (구간 [0,L]에서 최소화)
	if (!bSuccess)
	{
		const double phi = 0.5 * (3.0 - std::sqrt(5.0)); // 0.381966
		double a = 0.0, b = L;
		double x1 = a + (1.0 - phi) * (b - a);
		double x2 = a + phi * (b - a);
		double f1 = GetDist2BetweenStaPt(x1, q);
		double f2 = GetDist2BetweenStaPt(x2, q);

		for (int it = 0; it < 80; ++it) // 80회 반복
		{
			if (f1 > f2) 
			{
				a = x1;
				x1 = x2;
				f1 = f2;
				x2 = a + phi * (b - a);
				f2 = GetDist2BetweenStaPt(x2, q);
			}
			else 
			{
				b = x2;
				x2 = x1;
				f2 = f1;
				x1 = a + (1.0 - phi) * (b - a);
				f1 = GetDist2BetweenStaPt(x1, q);
			}

			if (std::fabs(b - a) < dTolS) break;
		}
		s = 0.5 * (a + b);
	}

	s = Clamp(s, 0.0, L);
	_Sta = s;

	return true;
}
