#include <iostream>
#include <iomanip>

#include "Structures.h"
#include "DataLoader.h"
#include "ShapeFunction.h"
#include "Jacobian.h"
#include "H.h"
#include "Hbc.h"
#include "P.h"
#include "t.h"
#include "C.h"

int main() {
	GlobalData gd;
	Grid grid;
	ElemUniv univ;

	//loadData("Test1_4_4.txt", gd, grid);
	loadData("Test2_4_4_MixGrid.txt", gd, grid);
	//loadData("Test3_31_31_kwadrat.txt", gd, grid);

	gd.npc = 2;
	int nPoints2D = gd.npc * gd.npc;

	for (int id : gd.BC) {
		grid.nodes[id - 1].BC = 1;
	}

	std::cout << std::fixed;

	std::cout << "Simulation Time: " << gd.SimulationTime << "\n";
	std::cout << "Number of Nodes*: " << gd.nN << "\n";
	std::cout << "Number of Elements*: " << gd.nE << "\n";
	
	std::cout << "BC* : ";
	for (int i : gd.BC) {
		std::cout << i << ", ";
	}

	// ========= Wypisanie elementow i id ich wierzcholkow ============

	std::cout << "\n\n*element" << endl;

	for (int j = 0; j < gd.nE; ++j) {
		std::cout << grid.elements[j].id << ", ";
		for (int i = 0; i < 4; ++i) {
			std::cout << grid.elements[j].nodeIds[i] << ", ";
		}
		std::cout << endl;
	}

	// ================= Wypisanie wspolrzednych wierzcholków ================

	std::cout << "\n*node" << endl;

	for (int j = 0; j < gd.nN; ++j) {
		std::cout << grid.nodes[j].id << ", " << std::setprecision(9) << grid.nodes[j].x << ", \t" << grid.nodes[j].y << " ";
		std::cout << endl;
	}

	std::cout << endl;

	// ======== Uruchomienie shape function do wyliczenia pochodnych funkcji ksztaltu po ksi i eta ===========

	ShapeFunction(gd.npc, univ);

	std::cout << "\n\ndN/dKsi\n";
	for (int i = 0; i < nPoints2D; ++i) {
		for (int j = 0; j < 4; ++j)
			std::cout << setw(10) << fixed << setprecision(6) << univ.dN_dKsi[i][j];
		std::cout << endl;
	}

	std::cout << "\n\ndN/dEta\n";
	for (int i = 0; i < nPoints2D; ++i) {
		for (int j = 0; j < 4; ++j)
			std::cout << setw(10) << fixed << setprecision(6) << univ.dN_dEta[i][j];
		std::cout << endl;

	}

	std::cout << "\n\nN\n";
	for (int i = 0; i < nPoints2D; ++i) {
		for (int j = 0; j < 4; ++j)
			std::cout << setw(10) << fixed << setprecision(6) << univ.N[i][j];
		std::cout << endl;
	}

	// =========== Liczymy Jakobiany ============

	computeJacobiansForAll(grid, univ, gd.npc);

	for (const Element& elem : grid.elements) {
		cout << "\nMacierze Jakobianu dla elementu " << elem.id << ":\n";
		for (size_t p = 0; p < elem.jacobians.size(); ++p) {
			const Jakobian& J = elem.jacobians[p];
			cout << "  PC" << p + 1 << ":\n";
			cout << "    [" << J.J[0][0] << "  " << J.J[0][1] << "]\n";
			cout << "    [" << J.J[1][0] << "  " << J.J[1][1] << "]\n";
			cout << "detJ = " << J.detJ << "\n\n";
		}
	}

	// ============= dN_dx dN_dy, H lokalne, Hbc, wektor P i C ===================

	computeHForAll(grid, univ, gd.npc, gd.Conductivity);
	computeHbcForAll(grid, grid, univ, gd.npc, gd.Alfa);
	computePForAll(grid, grid, univ, gd.npc, gd.Alfa, gd.Tot);
	computeCForAll(grid, univ, gd.npc, gd.Density, gd.SpecificHeat);

	for (const Element e : grid.elements) {
		std::cout << "\n=====Element " << e.id << "=====\n";
		std::cout << "\ndN/dx\n";
		for (int i = 0; i < nPoints2D; ++i) {
			for (int j = 0; j < 4; ++j)
				std::cout << setw(10) << fixed << setprecision(6) << e.dN_dx[i][j];
			std::cout << endl;
		}

		std::cout << "\ndN/dy\n";
		for (int i = 0; i < nPoints2D; ++i) {
			for (int j = 0; j < 4; ++j)
				std::cout << setw(10) << fixed << setprecision(6) << e.dN_dy[i][j];
			std::cout << endl;
		}

		std::cout << "\n Macierz H :\n";
		for (int a = 0; a < 4; ++a) {
			for (int b = 0; b < 4; ++b)
				std::cout << setw(12) << setprecision(6) << fixed << e.H[a][b];
			std::cout << endl;
		}

		std::cout << "\n Macierz Hbc :\n";
		for (int a = 0; a < 4; ++a) {
			for (int b = 0; b < 4; ++b)
				std::cout << setw(12) << setprecision(6) << fixed << e.Hbc[a][b];
			std::cout << endl;
		}

		std::cout << "\n Wektor P :\n";
		for (int a = 0; a < 4; ++a) {
				std::cout << setw(12) << setprecision(6) << fixed << e.P[a] << endl;
		}

		std::cout << "\n Macierz C :\n";
		for (int a = 0; a < 4; ++a) {
			for (int b = 0; b < 4; ++b)
				std::cout << setw(12) << setprecision(6) << fixed << e.C[a][b];
			std::cout << endl;
		}
	}

	// ==================== Macierz HG i PG i CG ======================
	Equation eq;
	eq.resetHG(gd.nN);
	eq.resetPG(gd.nN);
	eq.resetCG(gd.nN);

	for (const Element e : grid.elements) {
		addLocalToHG(eq, e.nodeIds, e.H, e.Hbc);
		addLocalToPG(eq, e.nodeIds, e.P);
		addLocalToCG(eq, e.nodeIds, e.C);
	}


	std::cout << "\n--- HG = H + HBC ---" << endl;

	for (int i = 0; i < gd.nN; ++i) {
		for (int j = 0; j < gd.nN; ++j) {
			std::cout << std::setw(7) << std::fixed << std::setprecision(3) << eq.HG[i][j];
		}
		std::cout << endl;
	}

	std::cout << "\n--- PG ---" << endl;

	for (int i = 0; i < gd.nN; ++i) {
		std::cout << std::setw(8) << std::fixed << std::setprecision(1) << eq.PG[i];
	}
	std::cout << endl;

	std::cout << "\n--- CG ---" << endl;

	for (int i = 0; i < gd.nN; ++i) {
		for (int j = 0; j < gd.nN; ++j) {
			std::cout << std::setw(9) << std::fixed << std::setprecision(3) << eq.CG[i][j];
		}
		std::cout << endl;
	}

	// ================== Rozwiazanie ukladu ==================
	std::cout << "\n\n========== Rozwiazanie ukladu ==========\n\n";

	eq.resetT(gd.nN);
	SolveT(eq, gd.nN, gd.InitialTemp, gd.SimulationStepTime, gd.SimulationTime, gd.Tot);

	std::cout << endl;
	for (int i = 0; i < gd.nN; ++i) {
		std::cout << "t[" << i + 1 << "] = " << std::fixed
			<< std::setprecision(6) << eq.t[i] << endl;
	}

	return 0;

}