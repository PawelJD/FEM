#pragma once
#include <iostream>
#include "Structures.h"

GaussLegendreData Data(int n);
double f1(double x);
double f2(double x, double y);
double gaussLegendre1D(double (*f)(double), int n, double a, double b);
double gaussLegendre2D(double (*f)(double, double), int n, double a, double b, double c, double d);
