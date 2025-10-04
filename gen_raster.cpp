#include <netcdf>
#include <iostream>
#include <vector>
#include <iostream>
#include <fstream>

using namespace std;
using namespace netCDF;
using namespace netCDF::exceptions;

int main(int argc, char* argv[]) {
    try {
        if (argc != 3) {
            cerr << "ERROR! USE " << argv[0] << " <source name> <target name> " << endl;
		    return -1;
        }

        string sourceName = argv[1]; // archivo de origen (texto plano)
        string targetName = argv[2]; // archivo de destino (.nc)

        // Abrir el archvio de origen
        ifstream in(sourceName);
        if (!in) {
            cerr << "Error abriendo el archivo de origen\n";
            return 1;
        }

        // ===================  Leer el archivo de origen  ===================

        // Se espera que los primeros dos numeros indiquen el tamaño del raster; filas x columnas
        int rows, columns;
        in >> rows >> columns;

        vector<double> raster(rows * columns); // raster usando indexación lineal

        // Rellenar el raster
        for (int i = 0; i < rows; i++)
        {
            for (int j = 0; j < columns; j++)
            {
                double value;
                in >> value; // lee el siguiente número en el archivo
                raster[i * columns + j] = value;
            }
        }

        // ===================  Escribir el archivo de destino  ===================


        NcFile newFile(targetName, NcFile::replace);

        //Definir dimensiones
        NcDim latDim = newFile.addDim("lat", rows);
        NcDim lonDim = newFile.addDim("lon", columns);

        // Crear el raster como variabe netcdf
        NcVar ncRaster = newFile.addVar("raster", ncDouble, {latDim, lonDim});

        // Escribirlos al archivo
        ncRaster.putVar(raster.data());

        cout << "Archivo " << targetName << " creado." << endl;

    } catch (NcException& e) {
        cerr << "Error: " << e.what() << endl;
        return 1;
    }
}
