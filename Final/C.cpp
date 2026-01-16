#include "C.h"
#include "GaussIntegration.h" 

using namespace std;

void computeCForElement(Element& elem, const Grid&, const ElemUniv& univ, int npc, double ro, double Cp)
{
    GaussLegendreData data = Data(npc);
    const int nPC = npc * npc;
 
    for (int a = 0; a < 4; ++a)
        for (int b = 0; b < 4; ++b)
            elem.C[a][b] = 0.0;

    // dla ka¿dego punktu ca³kowania p
    for (int p = 0; p < nPC; ++p) {
        const Jakobian& J = elem.jacobians[p];

        // waga 2D dla tego punktu: w_i * w_j
        int i = p / npc;     // indeks po ksi
        int j = p % npc;     // indeks po eta
        double w = data.w[i] * data.w[j];

        double detJ = std::abs(J.detJ);

        for (int a = 0; a < 4; ++a)
            for (int b = 0; b < 4; ++b)
                elem.C[a][b] += ro * Cp * univ.N[p][a] * univ.N[p][b] * detJ * w;
    }
}

void computeCForAll(Grid& grid, const ElemUniv& univ, int npc, double ro, double Cp)
{
    for (Element& elem : grid.elements)
        computeCForElement(elem, grid, univ, npc, ro, Cp);
}

void addLocalToCG(Equation& eq, const int nodeIds[4], const double Cloc[4][4]) {
    int I[4];
    for (int a = 0; a < 4; ++a) {
        I[a] = nodeIds[a] - 1;
    }
    for (int a = 0; a < 4; ++a)
        for (int b = 0; b < 4; ++b)
            eq.CG[I[a]][I[b]] += Cloc[a][b];
}
