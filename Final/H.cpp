#include "H.h"
#include "Hbc.h"
#include "GaussIntegration.h" 

using namespace std;

void computeHForElement(Element& elem, const Grid&, const ElemUniv& univ, int npc, double k) {
	// wagi Gaussa dla 1D, póŸniej waga 2D jako iloczyn
	GaussLegendreData data = Data(npc);
	const int nPC = npc * npc;

	elem.dN_dx.resize(nPC);
	elem.dN_dy.resize(nPC);

	// zerowanko macierzy H
	for (int a = 0; a < 4; ++a) {
		for (int b = 0; b < 4; ++b) {
			elem.H[a][b] = 0.0;
		}
	}

	// dla kazdego punktu calkowania
	for (int p = 0; p < nPC; ++p) {
		const Jakobian& J = elem.jacobians[p];

		for (int a = 0; a < 4; ++a) {
			double dNdKsi = univ.dN_dKsi[p][a];
			double dNdEta = univ.dN_dEta[p][a];
			elem.dN_dx[p][a] = J.invJ[0][0] * dNdKsi + J.invJ[1][0] * dNdEta;
			elem.dN_dy[p][a] = J.invJ[0][1] * dNdKsi + J.invJ[1][1] * dNdEta;
		}

		// waga 2D dla tego punktu: w_i * w_j
		int j = p / npc;     // indeks po ksi
		int i = p % npc;     // indeks po eta
		double w = data.w[i] * data.w[j];

		double detJ = std::abs(J.detJ);

		// sk³adka do H
		for (int a = 0; a < 4; ++a)
			for (int b = 0; b < 4; ++b)
				elem.H[a][b] += k * (elem.dN_dx[p][a] * elem.dN_dx[p][b] + elem.dN_dy[p][a] * elem.dN_dy[p][b]) * detJ * w;
	}


}

void computeHForAll(Grid& grid, const ElemUniv& univ, int npc, double k)
{
	for (Element& elem : grid.elements)
		computeHForElement(elem, grid, univ, npc, k);
}

void addLocalToHG(Equation& eq, const int nodeIds[4], const double Hloc[4][4], const double Hbcloc[4][4]) {
	int I[4];
	for (int a = 0; a < 4; ++a) {
		I[a] = nodeIds[a] - 1;
	}
	for (int a = 0; a < 4; ++a)
		for (int b = 0; b < 4; ++b)
			eq.HG[I[a]][I[b]] += (Hloc[a][b] + Hbcloc[a][b]);
}