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
    // 징후성: C(-x) = -C(x), S(-x) = -S(x)
    const double sign = (x < 0.0) ? -1.0 : 1.0;
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
    // u = s / (A * sqrt(pi))
    const double u = s / (A * std::sqrt(M_PI));

    double C, S;
    fresnelCS(u, C, S); // 표준 Fresnel C,S (π/2 t^2 정의)

    // 로컬 좌표 증가량 (시작 접선각 기준)
    const double scale = A * std::sqrt(M_PI);
    Vector2d dL{ scale * C, scale * S };

    // 시작 접선각 th0 만큼 회전 후, 시작점 p0에 더함
    Vector2d dG = dL.rotate(th0);
    return { p0.x + dG.x, p0.y + dG.y };
}

bool Clothoid::SetSByPoint(const Point2d& q)
{
    // ===== 설정값/유틸 =====
    auto clamp = [](double x, double a, double b) 
    {
        return x < a ? a : (x > b ? b : x);
    };

    // 최대 허용 길이 (Arc와 유사하게 _MaxSta를 쓴다고 가정)
    // 프로젝트에 맞게 교체 필요: 없으면 멤버/게터에서 가져오세요.
    const double L = /* _MaxSta */ this->_Sta; // <-- 실제 최대 길이로 교체 권장
    if (L <= 0.0) { _Sta = 0.0; return false; }

    // 현재 곡선의 기하

    auto Tvec = [&](double s) -> Vector2d {
        double th = theta(s);
        return { std::cos(th), std::sin(th) }; // unit tangent
        };
    auto Nvec = [&](double s) -> Vector2d {
        double th = theta(s);
        return { -std::sin(th), std::cos(th) }; // unit normal (좌수/우수 주의)
        };
    auto kappa = [&](double s) { return s / (A * A); }; // curvature(s)

    auto pos = [&](double s) -> Point2d { return this->position(s); };

    auto sqDist = [&](double s) {
        Point2d p = pos(s);
        double dx = p.x - q.x;
        double dy = p.y - q.y;
        return dx * dx + dy * dy;
        };

    // f(s) = ||p(s)-q||^2 의 도함수/이도함수
    auto derivs = [&](double s, double& fp, double& fpp) {
        Point2d p = pos(s);
        Vector2d r{ p.x - q.x, p.y - q.y };
        Vector2d T = Tvec(s);
        Vector2d N = Nvec(s);
        double k = kappa(s);
        // f'(s) = 2 r·T
        fp = 2.0 * (r.x * T.x + r.y * T.y);
        // f''(s) = 2( |T|^2 + r·(k N) ) = 2(1 + k (r·N))  (T는 단위벡터)
        fpp = 2.0 * (1.0 + k * (r.x * N.x + r.y * N.y));
        };

    // ===== 초기값 (시작점 접선으로 직선 투영) =====
    Vector2d T0{ std::cos(th0), std::sin(th0) };
    Vector2d r0{ q.x - p0.x, q.y - p0.y };
    double s = clamp(r0.x * T0.x + r0.y * T0.y, 0.0, L);

    // ===== 뉴턴법 =====
    const int    kMaxNewton = 12;
    const double tol_s = 1e-10 * std::max(1.0, L);
    const double tol_grad = 1e-10;

    bool newton_ok = false;
    for (int it = 0; it < kMaxNewton; ++it) {
        double fp, fpp;
        derivs(s, fp, fpp);

        if (std::fabs(fp) < tol_grad) { // 거의 최솟점
            newton_ok = true;
            break;
        }
        if (fpp == 0.0) break; // 비정상

        double step = fp / fpp;
        double s_new = s - step;

        // 구간 밖이면 클램프
        if (s_new < 0.0 || s_new > L) {
            // 가드: 경계 쪽이 더 낫다면 경계로 이동
            double f_here = sqDist(s);
            double f0 = sqDist(0.0);
            double fL = sqDist(L);
            if (f0 < f_here || fL < f_here) {
                s = (f0 < fL) ? 0.0 : L;
                newton_ok = true;
                break;
            }
            // 아니면 소폭 줄여서 시도
            s_new = clamp(s_new, 0.0, L);
        }

        if (std::fabs(s_new - s) < tol_s) {
            s = s_new;
            newton_ok = true;
            break;
        }
        s = s_new;
    }

    // ===== 백업: 골든 섹션 (구간 [0,L]에서 최소화) =====
    if (!newton_ok) {
        const double phi = 0.5 * (3.0 - std::sqrt(5.0)); // ≈ 0.381966
        double a = 0.0, b = L;
        double x1 = a + (1.0 - phi) * (b - a);
        double x2 = a + phi * (b - a);
        double f1 = sqDist(x1);
        double f2 = sqDist(x2);

        for (int it = 0; it < 80; ++it) {
            if (f1 > f2) {
                a = x1;
                x1 = x2;
                f1 = f2;
                x2 = a + phi * (b - a);
                f2 = sqDist(x2);
            }
            else {
                b = x2;
                x2 = x1;
                f2 = f1;
                x1 = a + (1.0 - phi) * (b - a);
                f1 = sqDist(x1);
            }
            if (std::fabs(b - a) < tol_s) break;
        }
        s = 0.5 * (a + b);
    }

    // 최종 세팅 (경계 안정화)
    s = clamp(s, 0.0, L);
    _Sta = s;
    return true;
}

