#include <fstream>
#include <sstream>
#include <iostream>

#include "DataLoader.h"

void loadData(const std::string& filename, GlobalData& gd, Grid& grid) {
	std::ifstream file(filename);
	std::string line;

	while (getline(file, line)) {
		if (line.find("SimulationTime") != std::string::npos)
			gd.SimulationTime = stoi(line.substr(15));
		else if (line.find("SimulationStepTime") != std::string::npos)
			gd.SimulationStepTime = stoi(line.substr(19));
		else if (line.find("Conductivity") != std::string::npos)
			gd.Conductivity = stoi(line.substr(12));
		else if (line.find("Alfa") != std::string::npos)
			gd.Alfa = stoi(line.substr(5));
		else if (line.find("Tot") != std::string::npos)
			gd.Tot = stoi(line.substr(4));
		else if (line.find("InitialTemp") != std::string::npos)
			gd.InitialTemp = stoi(line.substr(12));
		else if (line.find("Density") != std::string::npos)
			gd.Density = stoi(line.substr(8));
		else if (line.find("SpecificHeat") != std::string::npos)
			gd.SpecificHeat = stoi(line.substr(13));
		else if (line.find("Nodes number") != std::string::npos)
			gd.nN = stoi(line.substr(13));
		else if (line.find("Elements number") != std::string::npos)
			gd.nE = stoi(line.substr(16));
		else if (line.find("*Node") != std::string::npos) {
			for (int i = 0; i < gd.nN; ++i) {
				getline(file, line);
				std::istringstream part(line);
				Node node;
				part >> node.id;
				part.ignore(1, ',');
				part >> node.x;
				part.ignore(1, ',');
				part >> node.y;
				grid.nodes.push_back(node);

			}
		}
		else if (line.find("*Element") != std::string::npos) {
			for (int i = 0; i < gd.nE; ++i) {
				getline(file, line);
				std::istringstream part(line);
				Element elem;
				part >> elem.id;
				part.ignore(1, ',');
				for (int j = 0; j < 4; ++j) {
					part >> elem.nodeIds[j];
					if (j < 3) part.ignore(1, ',');
				}
				grid.elements.push_back(elem);
			}
		}
		else if (line.find("*BC") != std::string::npos) {
			getline(file, line);
			std::istringstream part(line);
			int liczba;
			while (!part.eof()) {
				part >> liczba;
				part.ignore(1, ',');
				gd.BC.push_back(liczba);
			}

		}
	}
}