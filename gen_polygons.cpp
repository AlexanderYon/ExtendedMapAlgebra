#include "ogrsf_frmts.h"
#include <fstream>
#include <vector>
#include <iostream>

struct Point { double x, y; };

int main(int argc, char** argv) {
	
	if(argc != 3){
		std::cerr << "ERROR! USE " << argv[0] << " <txtPath> <shpPath> " << std::endl;
		return -1;
	}
	
    OGRRegisterAll();
	const char* txtPath = argv[1]; // ruta archivo fuente
    const char* shpPath = argv[2]; // ruta archivo objetivo

    GDALDriver* poDriver = GetGDALDriverManager()->GetDriverByName("ESRI Shapefile");
    GDALDataset* poDS = poDriver->Create(shpPath, 0, 0, 0, GDT_Unknown, NULL);

    OGRLayer* poLayer = poDS->CreateLayer("polygons", nullptr, wkbPolygon, NULL);

    // Crear campo ID
    OGRFieldDefn oField("ID", OFTInteger);
    poLayer->CreateField(&oField);

    std::ifstream file(txtPath);
    if (!file.is_open()) {
        std::cerr << "No se pudo abrir el archivo\n";
        return 1;
    }

    int nPoligonos;
    file >> nPoligonos;

    for (int i = 0; i < nPoligonos; ++i) {
        int id, n_v;
        file >> id >> n_v;

        std::vector<Point> vertices(n_v);
        for (int j = 0; j < n_v; ++j) {
            file >> vertices[j].x >> vertices[j].y;
        }

        // Crear polígono
        OGRPolygon polygon;
        OGRLinearRing ring;
        for (auto& v : vertices)
            ring.addPoint(v.x, v.y);
        ring.closeRings();  // cierra automáticamente el anillo
        polygon.addRing(&ring);

        // Crear feature
        OGRFeature* poFeature = OGRFeature::CreateFeature(poLayer->GetLayerDefn());
        poFeature->SetField("ID", id);
        poFeature->SetGeometry(&polygon);
        poLayer->CreateFeature(poFeature);
        OGRFeature::DestroyFeature(poFeature);
    }

    GDALClose(poDS);
    std::cout << "Shapefile creado correctamente con " << nPoligonos << " polígonos.\n";
    return 0;
}
