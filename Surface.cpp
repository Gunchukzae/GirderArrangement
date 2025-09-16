#include "Surface.h"
#include "Base.h"

Surface::Surface()
{
	Init();
}

Surface::~Surface()
{
	RemoveAll();
}

void Surface::Init()
{
	ptOrigin = Point3d();
	sx = 0.0;
	sy = 0.0;
	sxx = 0.0;
	syy = 0.0;
	sxy = 0.0;

	pPlanarCrv = nullptr;
	pLongi = nullptr;
	pCross = nullptr;
}

void Surface::RemoveAll()
{
	if (pPlanarCrv != nullptr)
	{
		delete pPlanarCrv;
		pPlanarCrv = nullptr;
	}

	if (pLongi != nullptr)
	{
		delete pLongi;
		pLongi = nullptr;
	}

	if (pCross != nullptr)
	{
		delete pCross;
		pCross = nullptr;
	}
}

void Surface::ComputeSecondDerivsFromMeta()
{
	const double tol = 1e-9;

	// 종단면만 사용 (필수)
	if (pLongi) 
	{
		sx = pLongi->slope(ptOrigin.x);
		sxx = pLongi->sxx();
		// (선택) 원점 고도 동기화
		// ptOrigin.z = pLongi->elevation(ptOrigin.x);
	}
	else 
	{
		// 종단면이 반드시 있어야 한다면, 안전하게 0으로 두거나 assert/log
		sx = 0.0; sxx = 0.0;
		// assert(false && "pLongi must be set");
	}

	// Transverse curvature syy
	syy = 0.0; // always Line

	// Twist sxy
	sxy = sy * sxx; // 종단 곡률 존재 시
}

double Surface::Elevation(double x, double y) const
{
	const double dx = x - ptOrigin.x;
	const double dy = y - ptOrigin.y;

	double z0 = ptOrigin.z;
	double dz_longi = 0.0;

	if (pLongi)
	{
		// 종단면 기준의 정확한 프로파일 차분 사용
		const double z_at_x = pLongi->elevation(x);
		const double z_at_x0 = pLongi->elevation(ptOrigin.x);

		// 원점 고도를 종단면 기준으로 동기화했다면(권장), ptOrigin.z == z_at_x0
		// 동기화를 안 했더라도 차분만 더하면 double count가 없음
		dz_longi = (z_at_x - z_at_x0);
		z0 = pLongi->elevation(ptOrigin.x); // 원점 고도를 종단 기준으로
		// 또는 z0 = ptOrigin.z; 로 두고 dz_longi만 사용해도 됨(일관성만 유지)
	}
	else
	{
		// 종단면이 없으면 현행 2차식 근사
		dz_longi = sx * dx + 0.5 * sxx * dx * dx;
	}

	// 횡단(Line)과 twist
	const double dz_cross = sy * dy + 0.5 * syy * dy * dy + sxy * dx * dy;

	return z0 + dz_longi + dz_cross;
}

Surface::Type Surface::GetType() const
{
	const bool zero_sxx = (std::abs(sxx) < tol);
	const bool zero_syy = (std::abs(syy) < tol);
	const bool zero_sxy = (std::abs(sxy) < tol);

	if (zero_sxx && zero_syy && zero_sxy) return Type::Plane;
	if (zero_sxy) return Type::CurvatureNoTwist;
	return Type::CurvatureWithTwist;
}
