#include "Base.h"
#include "Surface.h"

#include <iostream>    
#include <vector>

static const char* AlignName(Alignment::Type t)
{
	switch (t)
	{
	case Alignment::Type::Line:      return "Alignment: Line";
	case Alignment::Type::Arc:       return "Alignment: Arc";
	case Alignment::Type::Clothoid:  return "Alignment: Clothoid";
	case Alignment::Type::Parabola3D:return "Alignment: Parabola3D";
	}

	return "Alignment: ?";
}
static const char* LongiName(LongitudinalSection::Type t)
{
	switch (t)
	{
	case LongitudinalSection::Type::Line:     return "Longitudinal: Line";
	case LongitudinalSection::Type::Arc:      return "Longitudinal: Arc";
	case LongitudinalSection::Type::Parabola2D: return "Longitudinal: Parabola2D";
	}

	return "Longitudinal: ?";
}
static const char* SrfTypeName(Surface::Type t)
{
	switch (t)
	{
	case Surface::Type::Plane:             return "Surface::Plane";
	case Surface::Type::CurvatureNoTwist:  return "Surface::CurvatureNoTwist";
	case Surface::Type::CurvatureWithTwist:return "Surface::CurvatureWithTwist";
	}

	return "Surface::?";
}

int main() {
	// Common
	auto pPlanarCrv = std::make_unique<IPlanarCurve>();
	auto pCross = std::make_unique<CrossSection>();
	auto pLongi = std::make_unique<LongitudinalSection>();
	pCross->type = CrossSection::Type::Line; // 항상 Line
	
	// Test
	Point3d origin(0.0, 0.0, 100.0);	// 원점
	const double sx_line = 1.0 / 40.0;  // 2.5% 종단
	const double sy_line = 1.0 / 20.0;  // 5%   횡단
	const Point3d query(2.0, 1.0, 0.0); // 시험점

	// Lambda
	auto run_case = [&](Alignment::Type alignType, LongitudinalSection::Type longiType,
		double Rv = 0.0, // Arc
		double L = 0.0, double G1 = 0.0, double G2 = 0.0)   // Parabola2D
	{
		//pPlanarCrv->type = alignType;

		// LongitudinalSection
		pLongi->type = longiType;
		pLongi->x0 = 0.0;
		pLongi->z0 = origin.z;

		switch (longiType) 
		{
		case LongitudinalSection::Type::Line:
			pLongi->sx0 = sx_line;     // 일정 경사
			break;
		case LongitudinalSection::Type::Arc:
			pLongi->sx0 = sx_line;     // 기준점 기울기
			pLongi->Rv = Rv;           // +: 오목, -: 볼록
			break;
		case LongitudinalSection::Type::Parabola2D:
			pLongi->L = L;
			pLongi->G1 = G1;          // 시작 기울기
			pLongi->G2 = G2;          // 끝 기울기
			// (선택) 일관성 위해 sx0=G1로 두고 싶다면: pLongi->sx0 = G1;
			break;
		}

		// Surface
		Surface srf;
		srf.ptOrigin = origin;
		srf.sy = sy_line;			// 횡단 기울기
		srf.pPlanarCrv = pPlanarCrv.get();  // 표기용
		srf.pCross = pCross.get();  // 항상 Line
		srf.pLongi = pLongi.get();  // 종단면

		srf.ComputeSecondDerivsFromMeta();

		const double z = srf.Elevation(query.x, query.y);
		const auto st = srf.GetType();

		std::cout << "---------------------------------------------\n";
		std::cout << AlignName(alignType) << " | " << LongiName(longiType) << " | Cross: Line\n";
		std::cout << "sx=" << srf.sx
			<< ", sy=" << srf.sy
			<< ", sxx=" << srf.sxx
			<< ", syy=" << srf.syy
			<< ", sxy=" << srf.sxy << "\n";
		std::cout << "Type = " << SrfTypeName(st) << "\n";
		std::cout << "Z(2,1) = " << z << "\n";
	};

	// 1) Alignment: Line + (Longi: Line/Arc/Parabola2D)
	run_case(Alignment::Type::Line, LongitudinalSection::Type::Line);
	run_case(Alignment::Type::Line, LongitudinalSection::Type::Arc, /*Rv=*/ 500.0);
	run_case(Alignment::Type::Line, LongitudinalSection::Type::Parabola2D, /*Rv=*/0.0, /*L=*/200.0, /*G1=*/sx_line, /*G2=*/-0.01);

	// 2) Alignment: Arc + (Longi: Line/Arc/Parabola2D)
	run_case(Alignment::Type::Arc, LongitudinalSection::Type::Line);
	run_case(Alignment::Type::Arc, LongitudinalSection::Type::Arc, 500.0);
	run_case(Alignment::Type::Arc, LongitudinalSection::Type::Parabola2D, 0.0, 200.0, sx_line, -0.01);

	// 3) Alignment: Clothoid + (Longi: Line/Arc/Parabola2D)
	run_case(Alignment::Type::Clothoid, LongitudinalSection::Type::Line);
	run_case(Alignment::Type::Clothoid, LongitudinalSection::Type::Arc, 500.0);
	run_case(Alignment::Type::Clothoid, LongitudinalSection::Type::Parabola2D, 0.0, 200.0, sx_line, -0.01);

	// 4) Alignment: Parabola3D + (Longi: Line/Arc/Parabola2D)
	run_case(Alignment::Type::Parabola3D, LongitudinalSection::Type::Line);
	run_case(Alignment::Type::Parabola3D, LongitudinalSection::Type::Arc, 500.0);
	run_case(Alignment::Type::Parabola3D, LongitudinalSection::Type::Parabola2D, 0.0, 200.0, sx_line, -0.01);

	return 0;
}
