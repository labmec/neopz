#include <iostream>
#include <fstream>
#include <string>
#include "JSON.hpp"

using json = nlohmann::json;
std::vector<double> get_tri_center(const std::vector<double> &pt0, const std::vector<double> &pt1, const std::vector<double> &pt2) {
	double x = (pt0[0] + pt1[0] + pt2[0]) / 3;
	double y = (pt0[1] + pt1[1] + pt2[1]) / 3;
	double z = (pt0[2] + pt1[2] + pt2[2]) / 3;
	return { x,y,z };
}

std::vector<double> get_mid_pt(const std::vector<double> &pt0, const std::vector<double> &pt1) {
	double x = (pt0[0] + pt1[0]) / 2;
	double y = (pt0[1] + pt1[1]) / 2;
	double z = (pt0[2] + pt1[2]) / 2;
	return { x,y,z };
}

void export_to_vtk(const std::vector<std::vector<double>> &coor, 
	const std::vector<std::vector<int>> &nodes, const std::vector<int> &sc, const std::vector<int> &mat) {
	std::ofstream myfile;
	myfile.open("rect.vtk");
	myfile << "# vtk DataFile Version 3.1" << std::endl;
	myfile << "Something" << std::endl;
	myfile << "ASCII" << std::endl;
	myfile << "DATASET POLYDATA" << std::endl;

	int numPts = coor.size();
	myfile << "POINTS\t" << numPts << "\tFLOAT" << std::endl;
	for (std::vector<double> pt : coor) {
		myfile << std::fixed << std::setprecision(5);
		pt[0] > -1e-6 ? myfile << "\t " << pt[0] : myfile << "\t" << pt[0];
		pt[1] > -1e-6 ? myfile << "\t\t " << pt[1] : myfile << "\t\t" << pt[1];
		pt[2] > -1e-6 ? myfile << "\t\t " << pt[2] : myfile << "\t\t" << pt[2];
		myfile << std::endl;
	}
	myfile << std::endl;

	int numCells = nodes.size();
	int numBd = sc.end() - std::find(sc.begin(), sc.end(), -1);
	myfile << "POLYGONS\t\t" << numCells-numBd << "\t\t" << (numCells- numBd )/5 * 21 << std::endl;
	int currentColor = 0;
	std::vector<int> cellColor;
	for (size_t i = 0; i < nodes.size(); i++) {
		if (mat[i] < 0) continue;
		std::vector<int> _faces = nodes[i];
		myfile << "\t\t" << _faces.size();
		for (int idx : _faces) {
			myfile << "\t\t" << idx;
		}
		myfile << std::endl;
		cellColor.push_back(mat[i]);
	}
	//for (size_t i = 0; i < nodes.size(); i++) {
	//	myfile << "\t\t" << 1;
	//	myfile << "\t\t" << sc[i];
	//	myfile << std::endl;
	//	cellColor.push_back(0);
	//}
	myfile << std::endl;
	myfile << std::endl;

	myfile << "CELL_DATA\t" << numCells - numBd << std::endl;
	myfile << "SCALARS\tMat\tint" << std::endl;
	myfile << "LOOKUP_TABLE default" << std::endl;
	for (int j : cellColor)
		j > 0 ? myfile << "\t\t" << j << std::endl: myfile <<"";
	myfile.close();

	// sc
	myfile.clear();
	myfile.open("rect_sc.vtk");
	myfile << "# vtk DataFile Version 3.1" << std::endl;
	myfile << "Something" << std::endl;
	myfile << "ASCII" << std::endl;
	myfile << "DATASET UNSTRUCTURED_GRID" << std::endl;
	
	numCells -= numBd;
	myfile << "POINTS\t" << numCells << "\tFLOAT" << std::endl;
	for (int idx : sc) {
		if (idx == -1) continue;
		std::vector<double> pt = coor[idx];
		myfile << std::fixed << std::setprecision(5);
		pt[0] > -1e-6 ? myfile << "\t " << pt[0] : myfile << "\t" << pt[0];
		pt[1] > -1e-6 ? myfile << "\t\t " << pt[1] : myfile << "\t\t" << pt[1];
		pt[2] > -1e-6 ? myfile << "\t\t " << pt[2] : myfile << "\t\t" << pt[2];
		myfile << std::endl;
	}
	myfile << std::endl;

	myfile << "CELLS\t\t" << numCells << "\t\t" << numCells *2 << std::endl;
	for (size_t i = 0; i < numCells; i++) {
		myfile << "\t\t" << 1;
		myfile << "\t\t" << i;
		myfile << std::endl;
	}
	myfile << std::endl;
	myfile << "CELL_TYPES\t" << numCells << std::endl;
	for (int i = 0; i < numCells; i++) myfile << "1" << std::endl;
	myfile << std::endl;
	myfile << "POINT_DATA " << numCells << std::endl;
	myfile << "SCALARS\tnode_id\tint"  << std::endl;
	myfile << "LOOKUP_TABLE\tdefault"  << std::endl;
	for (int i = 0; i < numCells; i++) myfile << i << std::endl;
	myfile.close();


	// bd
	myfile.clear();
	myfile.open("rect_bd.vtk");
	myfile << "# vtk DataFile Version 3.1" << std::endl;
	myfile << "Something" << std::endl;
	myfile << "ASCII" << std::endl;
	myfile << "DATASET UNSTRUCTURED_GRID" << std::endl;
	numCells -= numBd;
	myfile << "POINTS\t" << numPts << "\tFLOAT" << std::endl;
	for (std::vector<double> pt : coor) {
		myfile << std::fixed << std::setprecision(5);
		pt[0] > -1e-6 ? myfile << "\t " << pt[0] : myfile << "\t" << pt[0];
		pt[1] > -1e-6 ? myfile << "\t\t " << pt[1] : myfile << "\t\t" << pt[1];
		pt[2] > -1e-6 ? myfile << "\t\t " << pt[2] : myfile << "\t\t" << pt[2];
		myfile << std::endl;
	}
	myfile << std::endl;

	myfile << "CELLS\t\t" << numBd << "\t\t" << numBd * 3 << std::endl;
	for (size_t i = 0; i < numBd; i++) {
		myfile << "\t\t" << 2;
		myfile << "\t\t" << (*(nodes.end() - numBd + i))[0] <<"\t\t" << (*(nodes.end() - numBd + i))[1];
		myfile << std::endl;
	}
	myfile << std::endl;
	myfile << "CELL_TYPES\t" << numBd << std::endl;
	for (int i = 0; i < numBd; i++) myfile << "3" << std::endl;
	myfile << std::endl;
	myfile << "CELL_DATA " << numBd << std::endl;
	myfile << "SCALARS\tnode_id\tint" << std::endl;
	myfile << "LOOKUP_TABLE\tdefault" << std::endl;
	for (int i = 0; i < numBd; i++) myfile << mat[mat.size()-numBd+i]  << std::endl;
	myfile.close();
	std::cout << "exported to VTK file" << std::endl;
}
void export_to_json(const std::vector<std::vector<double>> &coor,
	const std::vector<std::vector<int>> &nodes, const std::vector<int> &sc, const std::vector<int> &mat) {
	json out;
	out["coor"] = coor;
	std::vector<json> elems;
	for (size_t i = 0; i < sc.size();i++) {
		json ele;
		ele["sc"] = sc[i];
		ele["nodes"] = nodes[i];
		ele["mat"] = mat[i];
		elems.push_back(ele);
	}
	out["elem"] = elems;
	std::ofstream myfile;
	myfile.open("rect.json");
	myfile << out.dump();
	myfile.close();
}
void rect_mesh(int numEleVer, double vert_domainsize) {
	double eleSize = vert_domainsize/numEleVer;
	int numEleHor = numEleVer*2;
	// generating grid points
	std::vector<std::vector<double>> coor;
	for (int i = 0; i <= numEleHor * 2; i++) {
		double xCoor = i*eleSize / 2;
		for (int j = 0; j <= numEleVer * 2; j++) {
			coor.push_back({ xCoor, j*eleSize / 2,0 });
		}
	}
	// generating connecivity
	std::vector<std::vector<int>> nodes;
	std::vector<int> sc;
	std::vector<int> mat;
	for (int row = 0; row < numEleHor; row++) {
		for (int col = 0; col < numEleVer; col++) {
			std::vector<int> ptsIdx;
			for (int k = 0; k < 3; k++) for (int l = 0; l < 3;l++) ptsIdx.push_back((row*2+l)*(2 * numEleVer + 1) + 2 * col + k);
			for (int il = 0; il < 5;il++) mat.push_back(1+(row+col)%2);
			nodes.push_back({ ptsIdx[1],ptsIdx[5],ptsIdx[7],ptsIdx[3] }); sc.push_back(ptsIdx[4]);
			nodes.push_back({ ptsIdx[0],ptsIdx[1],ptsIdx[3] }); sc.push_back(ptsIdx[0]);
			nodes.push_back({ ptsIdx[1],ptsIdx[2],ptsIdx[5] }); sc.push_back(ptsIdx[2]);
			nodes.push_back({ ptsIdx[5],ptsIdx[8],ptsIdx[7] }); sc.push_back(ptsIdx[8]);
			nodes.push_back({ ptsIdx[3],ptsIdx[7],ptsIdx[6] }); sc.push_back(ptsIdx[6]);
		}
	}
	
	// adjusting scaling center
	for (int row = 0; row < numEleVer; row++) {
		for (int col = 0; col < numEleHor; col++) {
			std::vector<int> ptsIdx;
			for (int k = 0; k < 3; k++) for (int l = 0; l < 3; l++) ptsIdx.push_back((col * 2 + l)*(2 * numEleVer + 1) + 2 * row + k);
			std::vector<double> center;
			int cellIdx;
			cellIdx = col*numEleVer + row;
			if (row == 0 && col == 0) {
				center = get_tri_center(coor[ptsIdx[0]], coor[ptsIdx[1]], coor[ptsIdx[3]]);
				coor.push_back(center);
				sc[cellIdx * 5 + 1] = coor.size() - 1;

				center = get_mid_pt(coor[ptsIdx[7]], coor[ptsIdx[6]]);
				coor.push_back(center);
				sc[cellIdx * 5 + 4] = coor.size() - 1;
				sc[(cellIdx + 1) * 5 + 1] = coor.size() - 1;

				center = get_mid_pt(coor[ptsIdx[5]], coor[ptsIdx[2]]);
				coor.push_back(center);
				sc[cellIdx * 5 + 2] = coor.size() - 1;
				sc[(cellIdx + numEleVer) * 5 + 1] = coor.size() - 1;

			}
			else if (row == 0 && col == numEleHor-1) {
				center = get_tri_center(coor[ptsIdx[1]], coor[ptsIdx[2]], coor[ptsIdx[5]]);
				coor.push_back(center);
				sc[cellIdx * 5 + 2] = coor.size() - 1;

				center = get_mid_pt(coor[ptsIdx[7]], coor[ptsIdx[8]]);
				coor.push_back(center);
				sc[cellIdx * 5 + 3] = coor.size() - 1;
				sc[(cellIdx + 1) * 5 + 2] = coor.size() - 1;
			}
			else if (row == numEleVer -1 && col == 0) {
				center = get_tri_center(coor[ptsIdx[3]], coor[ptsIdx[6]], coor[ptsIdx[7]]);
				coor.push_back(center);
				sc[cellIdx * 5 + 4] = coor.size() - 1;

				center = get_mid_pt(coor[ptsIdx[5]], coor[ptsIdx[8]]);
				coor.push_back(center);
				sc[cellIdx * 5 + 3] = coor.size() - 1;
				sc[(cellIdx + numEleVer) * 5 + 4] = coor.size() - 1;
			}
			else if (row == numEleVer - 1 && col == numEleHor - 1) {
				center = get_tri_center(coor[ptsIdx[7]], coor[ptsIdx[8]], coor[ptsIdx[5]]);
				coor.push_back(center);
				sc[cellIdx * 5 + 3] = coor.size() - 1;
			}
			else if (row == 0) {
				center = get_mid_pt(coor[ptsIdx[5]], coor[ptsIdx[2]]);
				coor.push_back(center);
				sc[cellIdx * 5 + 2] = coor.size() - 1;
				sc[(cellIdx + numEleVer) * 5 + 1] = coor.size() - 1;
			}
			else if (col == 0) {
				center = get_mid_pt(coor[ptsIdx[7]], coor[ptsIdx[6]]);
				coor.push_back(center);
				sc[cellIdx * 5 + 4] = coor.size() - 1;
				sc[(cellIdx + 1) * 5 + 1] = coor.size() - 1;
			}
			else if (row == numEleVer - 1) {
				center = get_mid_pt(coor[ptsIdx[5]], coor[ptsIdx[8]]);
				coor.push_back(center);
				sc[cellIdx * 5 + 3] = coor.size() - 1;
				sc[(cellIdx + numEleVer) * 5 + 4] = coor.size() - 1;
			}
			else if (col == numEleHor - 1) {
				center = get_mid_pt(coor[ptsIdx[7]], coor[ptsIdx[8]]);
				coor.push_back(center);
				sc[cellIdx * 5 + 3] = coor.size() - 1;
				sc[(cellIdx + 1) * 5 + 2] = coor.size() - 1;
			}
			else
				continue;
		}
	}
	// add boundary
	for (int i = 0; i < numEleHor; i++) {
		for (int j = 0; j < 4;j++) sc.push_back(-1);
		mat.push_back(-1); mat.push_back(-1);
		std::vector<int> ptsIdx;
		for (int k = 0; k < 3; k++) for (int l = 0; l < 3; l++) ptsIdx.push_back((i * 2 + l)*(2 * numEleVer + 1) + k);
		nodes.push_back({ ptsIdx[0],ptsIdx[1] });
		nodes.push_back({ ptsIdx[1],ptsIdx[2] });
		mat.push_back(-3); mat.push_back(-3);
		ptsIdx.clear();
		for (int k = 0; k < 3; k++) for (int l = 0; l < 3; l++) ptsIdx.push_back((i * 2 + l)*(2 * numEleVer + 1) + 2*(numEleVer-1) + k);
		nodes.push_back({ ptsIdx[6],ptsIdx[7] });
		nodes.push_back({ ptsIdx[7],ptsIdx[8] });
	}
	for (int i = 0; i < numEleVer; i++) {
		for (int j = 0; j < 4; j++) sc.push_back(-1);
		mat.push_back(-4); mat.push_back(-4);
		std::vector<int> ptsIdx;
		for (int k = 0; k < 3; k++) for (int l = 0; l < 3; l++) ptsIdx.push_back(l*(2 * numEleVer + 1) + 2 * i + k);
		nodes.push_back({ ptsIdx[0],ptsIdx[3] });
		nodes.push_back({ ptsIdx[3],ptsIdx[6] });
		mat.push_back(-2); mat.push_back(-2);
		ptsIdx.clear();
		for (int k = 0; k < 3; k++) for (int l = 0; l < 3; l++) ptsIdx.push_back(((numEleHor-1) * 2 + l)*(2 * numEleVer + 1) + 2 * i + k);
		nodes.push_back({ ptsIdx[2],ptsIdx[5] });
		nodes.push_back({ ptsIdx[5],ptsIdx[8] });
	}
    // putting a point element to hold rigid body modes
    {
        mat.push_back(-5);
        nodes.push_back({0});
        sc.push_back(-1);
        mat.push_back(-6);
        nodes.push_back({numEleVer * 2});
        sc.push_back(-1);
    }

	export_to_vtk(coor,nodes,sc,mat);
	export_to_json(coor, nodes, sc, mat);
}

void tutorial() {	// Explain how to read data from a given json file.


	// read in json file
	std::ifstream myfile("test.json");
	json j;
	myfile >> j;

	// Part coor	
	std::vector<std::vector<double>> coor = j["coor"]; // "coor" are in 2d vector
	int nnodesTotal = coor.size();
	std::cout << "Coordinates ( " << nnodesTotal << " in total )" << std::endl;
	for (int i = 0; i < nnodesTotal; i++) {
		std::cout << "\t[";
		for (int j = 0; j < 2; j++)
			std::cout << coor[i][j] << " ,";
		std::cout << coor[i][2] << "]" << std::endl;
	}





	// Part elem
	std::vector<json> elem = j["elem"]; // "elem" are in 1d vector with json object
	int nElem = elem.size();
	std::cout << "Elements ( " << nElem << " in total)" << std::endl;
	for (int i = 0; i < nElem; i++) {
		std::cout << "\t Elem" << i << ":" << std::endl;
		// sc idx
		int sc = elem[i]["sc"]; //sc idx
		std::cout << "\t\t sc_idx: " << sc << std::endl;
		// node idx
		std::vector<int> nodes = elem[i]["nodes"]; // nodes list
		int nnodes = nodes.size();
		std::cout << "\t\t node_idx (" << nnodes << " in total ) :";
		for (int k = 0; k < nnodes - 1; k++)
			std::cout << nodes[k] << ",";
		std::cout << nodes[nnodes - 1] << std::endl;
		// mat idx
		int mat = elem[i]["mat"];
		std::cout << "\t\t mat_idx: " << mat << std::endl;
	}
}

int mainJSon() {
	//tutorial();
    int numver_el = 5;
    double domain_vertsize = 25.;
	rect_mesh(numver_el,domain_vertsize);
    return 0;
}

