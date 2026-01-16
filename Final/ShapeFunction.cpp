#include <array>

#include "ShapeFunction.h"
#include "GaussIntegration.h"

array<double, 4> Nvals(double ksi, double eta) {
	array<double, 4> N;
	N[0] = 0.25 * (1 - ksi) * (1 - eta);
	N[1] = 0.25 * (1 + ksi) * (1 - eta);
	N[2] = 0.25 * (1 + ksi) * (1 + eta);
	N[3] = 0.25 * (1 - ksi) * (1 + eta);
	return N;
}

void ShapeFunction(int npc, ElemUniv& univ) {
	GaussLegendreData data = Data(npc);
	vector<pair<double, double>> pc;
	int nPoints = npc * npc;

	univ.N.resize(nPoints);

	univ.dN_dKsi.resize(nPoints);
	univ.dN_dEta.resize(nPoints);

	// generowanie punktow calkowania
	for (int i = 0; i < npc; ++i) {
		for (int j = 0; j < npc; ++j) {
			pc.push_back({ data.x[i], data.x[j] });
		}
	}

	// pochodne funkcji ksztaltu w kazdym punkcie calkowania
	for (int i = 0; i < nPoints; ++i) {
		double ksi = pc[i].first;
		double eta = pc[i].second;

		// funkcja ksztaltu
		univ.N[i] = Nvals(ksi, eta);

		//pochodne wzgledem ksi
		univ.dN_dKsi[i][0] = -0.25 * (1 - eta);
		univ.dN_dKsi[i][1] = 0.25 * (1 - eta);
		univ.dN_dKsi[i][2] = 0.25 * (1 + eta);
		univ.dN_dKsi[i][3] = -0.25 * (1 + eta);

		//pochodne wzgledem eta
		univ.dN_dEta[i][0] = -0.25 * (1 - ksi);
		univ.dN_dEta[i][1] = -0.25 * (1 + ksi);
		univ.dN_dEta[i][2] = 0.25 * (1 + ksi);
		univ.dN_dEta[i][3] = 0.25 * (1 - ksi);
	}

	for (int e = 0; e < 4; ++e)
		univ.surfaces[e].N.resize(npc);

	for (int k = 0; k < npc; ++k) {
		double s = data.x[k];   // punkt Gaussa w 1D

		// krawêdŸ 0: dolna (1-2): eta = -1, ksi = s
		univ.surfaces[0].N[k] = Nvals(s, -1.0);

		// krawêdŸ 1: prawa (2-3): ksi = 1, eta = s
		univ.surfaces[1].N[k] = Nvals(1.0, s);

		// krawêdŸ 2: górna (3-4): eta = 1, ksi = s
		univ.surfaces[2].N[k] = Nvals(s, 1.0);

		// krawêdŸ 3: lewa (4-1): ksi = -1, eta = s
		univ.surfaces[3].N[k] = Nvals(-1.0, s);
	}

}