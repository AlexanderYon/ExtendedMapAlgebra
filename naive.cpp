#include <netcdf>
#include <iostream>
#include <vector>
#include <ogrsf_frmts.h>
#include <chrono>
#include <cstdio>
#include <iomanip>
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/polygon.hpp>

namespace bg = boost::geometry;
using BoostPoint = bg::model::d2::point_xy<double>;
using BoostPolygon = bg::model::polygon<BoostPoint>;
using namespace std;
using namespace netCDF;
using namespace netCDF::exceptions;
using namespace std::chrono;

/*

            HELPER FUNCTIONS

*/

/**
 * Get the MBR of a given Polygon
 */
std::tuple<int,int,int,int> getMBR(OGRPolygon* poly) {

    // Initial values
    double xmin = std::numeric_limits<double>::max();
    double xmax = std::numeric_limits<double>::lowest();
    double ymin = std::numeric_limits<double>::max();
    double ymax = std::numeric_limits<double>::lowest();

    // Iterate outer ring
    OGRLinearRing* ring = poly->getExteriorRing();
    int nPoints = ring->getNumPoints();

    for (int i = 0; i < nPoints; i++) {
        double x = ring->getX(i);
        double y = ring->getY(i);

        if (x < xmin) xmin = x;
        if (x > xmax) xmax = x;
        if (y < ymin) ymin = y;
        if (y > ymax) ymax = y;
    }

    // Raster origin is (0,0) and 1x1 resolution
    int ixmin = static_cast<int>(std::floor(xmin));
    int ixmax = static_cast<int>(std::ceil(xmax) - 1);
    int iymin = static_cast<int>(std::floor(ymin));
    int iymax = static_cast<int>(std::ceil(ymax) - 1);

    return {ixmin, ixmax, iymin, iymax};
}

// Convert the cell (i,j) of the raster to a Polygon
OGRPolygon fromCellToPolygon(int i, int j) {
    OGRPolygon cellPoly;
    OGRLinearRing ring;

    // Cell corners
    // The indices i and j are interchanged due to how geographic coordinates should be read.
    ring.addPoint(j, i);       
    ring.addPoint(j+1, i);     
    ring.addPoint(j+1, i+1);   
    ring.addPoint(j, i+1);     
    ring.closeRings();

    cellPoly.addRing(&ring);
    return cellPoly;
}


/**
 * Convert a OGRPolygon to a BoostPolygon
 * 
 */
BoostPolygon convertOGRPolygonToBoost(OGRPolygon* ogrPoly) {
    BoostPolygon boostPoly;

    // Anillo exterior
    OGRLinearRing* exteriorRing = ogrPoly->getExteriorRing();
    for (int i = 0; i < exteriorRing->getNumPoints(); ++i) {
        boostPoly.outer().emplace_back(exteriorRing->getX(i), exteriorRing->getY(i));
    }

    // Anillos interiores (agujeros)
    for (int i = 0; i < ogrPoly->getNumInteriorRings(); ++i) {
        OGRLinearRing* innerRing = ogrPoly->getInteriorRing(i);
        bg::model::ring<BoostPoint> boostInnerRing;
        for (int j = 0; j < innerRing->getNumPoints(); ++j) {
            boostInnerRing.emplace_back(innerRing->getX(j), innerRing->getY(j));
        }
        boostPoly.inners().push_back(boostInnerRing);
    }

    // Corrige y asegura que el polígono sea válido (cerrado, sentido correcto)
    bg::correct(boostPoly);
    return boostPoly;
}


void prinResultsTable(const std::map<long, double>& results) {
    cout << "RESULTS" << endl;
    cout << "| " << left << setw(20) << "Polygon ID" 
         << "| " << left << setw(20) << "Maximum" 
         << " |" << endl;
    for(const auto& pair : results) {
        cout << "| " << left << setw(20) << pair.first 
             << "| " << left << setw(20) << pair.second 
             << " |" << endl;
    }
}

/*

            MAIN

*/


int main(int argc, char* argv[]) {
    try {
        if (argc != 3) {
            cerr << "ERROR! USE " << argv[0] << " <nc_raster> <shapefile>" << endl;
		    return -1;
        }

        // Start Time
        using clock = high_resolution_clock;
        auto start_time = clock::now();

        // ================================     GET RASTER     ================================

        // Extract raster data
        NcFile dataFile(argv[1], NcFile::read);
        NcDim latDim = dataFile.getDim("lat"); // latitude
        NcDim lonDim = dataFile.getDim("lon"); // longitude
        size_t nLat  = latDim.getSize();       
        size_t nLon  = lonDim.getSize();

        // Extract raster
        NcVar raster = dataFile.getVar("raster");


        // ================================     GET POLYGONS     ================================

        // Open the .shp file to read
        const char* shpPath = argv[2];
        OGRRegisterAll();
        GDALDataset *poDS = (GDALDataset *)GDALOpenEx(shpPath, GDAL_OF_VECTOR, NULL, NULL, NULL);
        if (poDS == NULL) {
            cerr << "Could not open " << shpPath << endl;
            return 1;
        }

        // ================================     NAIVE ALGORITHM    ================================
        
        map<long, double> results; // results table: <PolygonID, Maximum>
        vector<long> excludedPolygonsIDs;
        int numberOfExcludedPolygons = 0;  // count how many polygons are excluded from the results table 

        // Get the first layer
        OGRLayer *layer = poDS->GetLayer(0);

        // Iterate through each feature (each polygon)
        OGRFeature *feature;
        layer->ResetReading();
        while ((feature = layer->GetNextFeature()) != NULL) {
            
            // Get the geometry of the current polygon
            OGRGeometry *geom = feature->GetGeometryRef();

            if (geom != NULL && wkbFlatten(geom->getGeometryType()) == wkbPolygon) {

                OGRPolygon *OGRPoly = (OGRPolygon *) geom;
                BoostPolygon poly = convertOGRPolygonToBoost(OGRPoly); // convert to BoostPolygon to use covered_by

                // Get the MBR of the current Polygon 
                auto [ixmin, ixmax, iymin, iymax] = getMBR(OGRPoly);

                // Define the raster subsection
                size_t start_row = iymin;
                size_t start_column = ixmin;
                size_t nrows = iymax - iymin + 1;
                size_t ncols = ixmax - ixmin + 1;
                vector<size_t> start = {start_row, start_column};   // cell (start_row, start_column)
                vector<size_t> count = {nrows, ncols};              // nrows rows, ncols columns

                // Get the raster subsection
                vector<float> subset(count[0] * count[1]);

                // Read subsection only
                raster.getVar(start, count, subset.data());
                double currentMax = numeric_limits<double>::lowest();

                // Cell-by-cell Comparisons
                for (size_t i = 0; i < count[0]; i++) {
                    for (size_t j = 0; j < count[1]; j++) {
                        OGRPolygon cellPolygon = fromCellToPolygon(i + iymin, j + ixmin);

                        BoostPolygon boostedCellPolygon = convertOGRPolygonToBoost(&cellPolygon);

                        if (bg::covered_by(boostedCellPolygon, poly)) {
                            double currentValue = subset[i * count[1] + j];
                            if (currentValue > currentMax) {
                                currentMax = currentValue;
                            }
                        }
                    }
                }

                // Save pair <PolygonID, Max> only if a maximum value was found
                if (currentMax > numeric_limits<double>::lowest()) {
                    results[feature->GetFID()] = currentMax;
                } else {
                    // Exclude the current Polygon
                    numberOfExcludedPolygons++;
                    excludedPolygonsIDs.push_back(feature->GetFID());
                }
            }
            // Destroy feature
            OGRFeature::DestroyFeature(feature);
        }
        GDALClose(poDS);

        // End time
        auto end_time = clock::now();

        // Calculate execution time
        auto execution_time = duration_cast<milliseconds>(end_time - start_time);

        // ================================     PRINT RESULTS    ================================
        
        prinResultsTable(results);
        cout << "Total execution time (including data loading): " << execution_time.count() << " ms "<< endl;
        cout << "Number of polygons excluded: " << numberOfExcludedPolygons << endl;
        cout << "Polygons excluded: ";

        for (const long& polygonID : excludedPolygonsIDs) {
            cout << polygonID << " ";
        }

        cout << endl;
        
        return 0;

    } catch (NcException& e) {
        cerr << "Error: " << e.what() << endl;
        return 1;
    }
}