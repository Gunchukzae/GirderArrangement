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

	// ���ܸ鸸 ��� (�ʼ�)
	if (pLongi) 
	{
		sx = pLongi->slope(ptOrigin.x);
		sxx = pLongi->sxx();
		// (����) ���� �� ����ȭ
		// ptOrigin.z = pLongi->elevation(ptOrigin.x);
	}
	else 
	{
		// ���ܸ��� �ݵ�� �־�� �Ѵٸ�, �����ϰ� 0���� �ΰų� assert/log
		sx = 0.0; sxx = 0.0;
		// assert(false && "pLongi must be set");
	}

	// Transverse curvature syy
	syy = 0.0; // always Line

	// Twist sxy
	sxy = sy * sxx; // ���� ��� ���� ��
}

double Surface::Elevation(double x, double y) const
{
	const double dx = x - ptOrigin.x;
	const double dy = y - ptOrigin.y;

	double z0 = ptOrigin.z;
	double dz_longi = 0.0;

	if (pLongi)
	{
		// ���ܸ� ������ ��Ȯ�� �������� ���� ���
		const double z_at_x = pLongi->elevation(x);
		const double z_at_x0 = pLongi->elevation(ptOrigin.x);

		// ���� ���� ���ܸ� �������� ����ȭ�ߴٸ�(����), ptOrigin.z == z_at_x0
		// ����ȭ�� �� �ߴ��� ���и� ���ϸ� double count�� ����
		dz_longi = (z_at_x - z_at_x0);
		z0 = pLongi->elevation(ptOrigin.x); // ���� ���� ���� ��������
		// �Ǵ� z0 = ptOrigin.z; �� �ΰ� dz_longi�� ����ص� ��(�ϰ����� ����)
	}
	else
	{
		// ���ܸ��� ������ ���� 2���� �ٻ�
		dz_longi = sx * dx + 0.5 * sxx * dx * dx;
	}

	// Ⱦ��(Line)�� twist
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
