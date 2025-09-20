#include "Clothoid.h"
#include "Base.h"
#include <corecrt_math_defines.h>
#include <limits>
#include <functional>
#include <algorithm>

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
    _MaxSta = 0.0;
    _Sta = 0.0;
}

void Clothoid::RemoveAll()
{

}

static double simpson(const std::function<double(double)>& f, double a, double b) 
{
    const double c = 0.5 * (a + b);

    return (b - a) * (f(a) + 4.0 * f(c) + f(b)) / 6.0;
}

static double adaptiveSimpson(const std::function<double(double)>& f, double a, double b, double eps, double whole, int depth, int maxDepth)
{
    const double c = 0.5 * (a + b);
    const double left = simpson(f, a, c);
    const double right = simpson(f, c, b);
    const double delta = left + right - whole;

    if (depth >= maxDepth || std::fabs(delta) < 15.0 * eps) 
        return left + right + delta / 15.0; // 리처드슨 보정

    return adaptiveSimpson(f, a, c, eps * 0.5, left, depth + 1, maxDepth) + adaptiveSimpson(f, c, b, eps * 0.5, right, depth + 1, maxDepth);
}

// Fresnel : C(x)=∫0^x cos(π/2 t^2) dt, S(x)=∫0^x sin(π/2 t^2) dt
static void fresnelCS(double x, double& C, double& S, double eps = 1e-10, int maxDepth = 20)
{
    const double sign = (x < 0.0) ? -1.0 : 1.0; // C,S는 홀함수
    const double ax = std::fabs(x);

    auto fc = [](double t) { return std::cos(0.5 * M_PI * t * t); };
    auto fs = [](double t) { return std::sin(0.5 * M_PI * t * t); };

    // 초기 심프슨 값
    const double Sc0 = simpson(fc, 0.0, ax);
    const double Ss0 = simpson(fs, 0.0, ax);

    const double Cpos = adaptiveSimpson(fc, 0.0, ax, eps, Sc0, 0, maxDepth);
    const double Spos = adaptiveSimpson(fs, 0.0, ax, eps, Ss0, 0, maxDepth);

    C = sign * Cpos;
    S = sign * Spos;
}

Point2d Clothoid::position(double s) const
{
    Point2d pt;

    const double u = s / (A * std::sqrt(M_PI)); // 치환
    double C, S;
    fresnelCS(u, C, S);

    const double scale = A * std::sqrt(M_PI); // 로컬 좌표 증가량
    Vector2d v0 = Vector2d(scale * C, scale * S);
    Vector2d v = v0.rotate(th0); // 시작 접선각만큼 회전
    pt = Point2d(p0.x + v.x, p0.y + v.y); // 시작점 이동

    return pt;
}

static double Clamp(double x, double a, double b)
{
    return x < a ? a : (x > b ? b : x);
}

Vector2d Clothoid::GetTanVecBySta(double s)
{
    double th = theta(s);

    return Vector2d(std::cos(th), std::sin(th));
}

Vector2d Clothoid::GetNormVecBySta(double s)
{
    double th = theta(s);

    return Vector2d(-std::sin(th), std::cos(th));
}

double Clothoid::GetDist2BetweenStaPt(const double s, const Point2d& q)
{
    Point2d p = position(s);
    double dx = p.x - q.x;
    double dy = p.y - q.y;

    return dx * dx + dy * dy;
}

void Clothoid::GetDistDeriveBetweenStaPt(double s, const Point2d& q, double& fp, double& fpp)
{
    Point2d p = position(s);
    Vector2d r = Vector2d(p.x - q.x, p.y - q.y); // r(s) = p(s) − q
    Vector2d T = GetTanVecBySta(s);
    Vector2d N = GetNormVecBySta(s);
    double k = curvature(s);
    
    // f(s) = r^2
    // f'(s) = 2 r·T
    fp = 2.0 * (r.x * T.x + r.y * T.y);

    // f''(s) = 2( |T|^2 + r·(k N) ) = 2(1 + k (r·N))
    fpp = 2.0 * (1.0 + k * (r.x * N.x + r.y * N.y));
}

bool Clothoid::SetStaByPoint(const Point2d& q)
{
    if (A != 0)
        _MaxSta = A * 5;
    
    const double L = _MaxSta;
    if (L <= 0.0) { _Sta = 0.0; return false; }

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

        if (fpp <= 0.0) break; // Error

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
