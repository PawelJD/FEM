#pragma once
#include <vector>
#include <array>

using namespace std;

struct Jakobian;

struct Node {
	int id;
	double x, y;
	int BC = 0;
};

struct  Element {
	int id;
	int nodeIds[4];
	vector<Jakobian> jacobians;
	vector<array<double, 4>> dN_dx;
	vector<array<double, 4>> dN_dy;
	double H[4][4] = {};
	double Hbc[4][4] = {};
	double P[4] = {};
	double C[4][4] = {};
};

struct GlobalData {
	int SimulationTime;
	int SimulationStepTime;
	int Conductivity;
	int Alfa;
	int Tot;
	double InitialTemp;
	int Density;
	int SpecificHeat;
	int nN; //liczba wezlow
	int nE; //liczba elementow
	vector<int> BC;
	int npc;

	GlobalData() : SimulationTime(0), SimulationStepTime(0), Conductivity(0), Alfa(0), Tot(0),
		InitialTemp(0.0), Density(0), SpecificHeat(0), nN(0), nE(0), npc(0) {}
};

struct Grid {
	vector<Node> nodes;
	vector<Element> elements;
};

struct ElemUniv {
	vector<array<double, 4>> N;

	vector<array<double, 4>> dN_dKsi;
	vector<array<double, 4>> dN_dEta;

	struct Surface {
		vector<array<double, 4>> N; // funkcje ksztaltu na powierzchni (dla kazdego punktu calkowania)
	};
	std::array<Surface, 4> surfaces; // 4 krawedzie
};

struct Jakobian {
	double J[2][2];
	double detJ;
	double invJ[2][2];
};

struct GaussLegendreData {
	const double* x;
	const double* w;
	int size;
};

struct Equation {
	vector<vector<double>> HG;
	vector<double> PG;
	vector<double> t;
	vector<vector<double>> CG;

	void resetHG(int n) {
		HG.assign(n, vector<double>(n, 0.0));
	}

	void resetPG(int n) {
		PG.assign(n, 0.0);
	}

	void resetT(int n) {
		t.assign(n, 0.0);
	}

	void resetCG(int n) {
		CG.assign(n, vector<double>(n, 0.0));
	}
};