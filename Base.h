#pragma once
#include <complex>

struct Point2d
{
	Point2d() : x(0.0), y(0.0) {}
	Point2d(double x_, double y_) : x(x_), y(y_) {}
	virtual ~Point2d() {}

	double x, y;
};

struct Point3d
{
	Point3d() : x(0.0), y(0.0), z(0.0) {}
	Point3d(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}
	virtual ~Point3d() {}

	double x, y, z;
};

struct Vector2d
{
	Vector2d() : x(0.0), y(0.0) {}
	Vector2d(double xx,	double yy) : x(xx), y(yy) {}
	virtual ~Vector2d() {}

	Vector2d rotate(double dth) const 
	{	
		double cos_t = std::cos(dth); 
		double sin_t = std::sin(dth);	
		return { x * cos_t - y * sin_t, x * sin_t + y * cos_t };
	}

	double x, y;
};
