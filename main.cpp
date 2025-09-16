#include <iostream>    
#include <vector>

#include "Base.h"
#include "Surface.h"
#include "Arc.h"

int main() {
	// Common
	auto pPlanarCrv = std::make_unique<Arc>();
	auto pLongi = std::make_unique<LongitudinalSection>();
	auto pCross = std::make_unique<CrossSection>();
	pCross->type = CrossSection::Type::Line; // �׻� Line
	
	// Test
	Point3d origin(0.0, 0.0, 100.0);	// ����
	const double sx_line = 1.0 / 40.0;  // 2.5% ����
	const double sy_line = 1.0 / 20.0;  // 5%   Ⱦ��
	const Point3d query(2.0, 1.0, 0.0); // ������

	return 0;
}
