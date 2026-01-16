#include "GaussIntegration.h"

using namespace std;

GaussLegendreData Data(int n) {
    static const double x2[] = { -(1.0 / sqrt(3.0)), 1.0 / sqrt(3.0) };
    static const double w2[] = { 1.0, 1.0 };

    static const double x3[] = { -(sqrt(3.0 / 5.0)), 0.0, sqrt(3.0 / 5.0) };
    static const double w3[] = { 5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0 };

    static const double x4[] = { -0.861136, -0.339981, 0.339981, 0.861136 };
    static const double w4[] = { 0.347855, 0.652145, 0.652145, 0.347855 };

    static const double x5[] = { -0.906180, -0.538469, 0.0, 0.538469, 0.906180 };
    static const double w5[] = { 0.236927, 0.478629, 0.568889, 0.478629, 0.236927 };

    switch (n) {
    case 2: return { x2, w2, 2 };
    case 3: return { x3, w3, 3 };
    case 4: return { x4, w4, 4 };
    case 5: return { x5, w5, 5 };
    default:
        cout << "Obslugiwane tylko 2 do 5 wezlow.\n";
        exit(1);
    }
}

// Funkcja 1: 5(x^2) + 3x + 6
double f1(double x) {
    return 5 * pow(x, 2) + 3 * x + 6;
}

// Funkcja 2: 5(x^2)(y^2) + 3xy + 6
double f2(double x, double y) {
    return 5 * pow(x, 2) * pow(y, 2) + 3 * x * y + 6;
}

// Kwadratura Gaussa-Legendre'a
double gaussLegendre1D(double (*f)(double), int n, double a, double b) {
    GaussLegendreData data = Data(n);
    double sum = 0.0;
    for (int i = 0; i < data.size; ++i) {
        double xi = data.x[i];
        double wi = data.w[i];
        double x = 0.5 * ((b - a) * xi + (b + a));
        sum += wi * f(x);
    }
    return 0.5 * (b - a) * sum;
}

double gaussLegendre2D(double (*f)(double, double), int n, double a, double b, double c, double d) {
    GaussLegendreData data = Data(n);
    double sum = 0.0;

    for (int i = 0; i < data.size; i++) {
        for (int j = 0; j < data.size; j++) {
            double xi = 0.5 * ((b - a) * data.x[i] + (b + a));
            double eta = 0.5 * ((d - c) * data.x[j] + (d + c));

            sum += data.w[i] * data.w[j] * f(xi, eta);
        }
    }

    return 0.25 * (b - a) * (d - c) * sum;
}