#include <netcdf>
#include <iostream>
#include <vector>
#include <ogrsf_frmts.h>
#include <chrono>
#include <cstdio>
#include <iomanip>
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/polygon.hpp>
#include <optional>

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
 * Convert an OGRPolygon to a BoostPolygon
 */
BoostPolygon convertOGRPolygonToBoost(OGRPolygon* ogrPoly) {
    BoostPolygon boostPoly;

    // Outer Ring
    OGRLinearRing* exteriorRing = ogrPoly->getExteriorRing();
    for (int i = 0; i < exteriorRing->getNumPoints(); ++i) {
        boostPoly.outer().emplace_back(exteriorRing->getX(i), exteriorRing->getY(i));
    }

    // Inner rings
    for (int i = 0; i < ogrPoly->getNumInteriorRings(); ++i) {
        OGRLinearRing* innerRing = ogrPoly->getInteriorRing(i);
        bg::model::ring<BoostPoint> boostInnerRing;
        for (int j = 0; j < innerRing->getNumPoints(); ++j) {
            boostInnerRing.emplace_back(innerRing->getX(j), innerRing->getY(j));
        }
        boostPoly.inners().push_back(boostInnerRing);
    }

    // Correct and make sure the Polygon is valid (closed, correct direction)
    bg::correct(boostPoly);
    return boostPoly;
}

/**
 * Print results table using this format:
RESULTS
| Polygon ID          | Maximum              |
| id1                 | max1                 |
| id2                 | max2                 |
| id3                 | max3                 |
|  .                  |  .                   |
|  .                  |  .                   |
|  .                  |  .                   |
| idn                 | maxn                 |
Total execution time (including data loading): <time> ms 
Number of polygons excluded: <number>
Polygons excluded: <id0> <id1> <id2> ... <idk>
 */
void printResultsTable(const std::map<long, double>& results) {
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

        // Define clock for time measuring
        using clock = high_resolution_clock;

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

        // Start time
        duration<double> total_netcdf_queries_time(0);      // measure NetCDF queries (for the Raster)
        duration<double> total_gdal_queries_time(0);        // measure GDAL queries (for Polygons)
        duration<double> total_extra_opreations_time(0);    // extra operations like conversions between objects
        duration<double> total_execution_time(0);
        int netcdfIOCount = 0;
        int gdalIOCount = 0;
        
        auto start_total_time = clock::now();
        while (true) {

            auto start_gdal_query_time = clock::now();
            feature = layer->GetNextFeature();
            auto end_gdal_query_time = clock::now();
            total_gdal_queries_time += (end_gdal_query_time - start_gdal_query_time);

            if (feature == NULL) {
                total_gdal_queries_time -= (end_gdal_query_time - start_gdal_query_time);
                break;
            }

            gdalIOCount++;
            
            // Get the geometry of the current polygon
            OGRGeometry *geom = feature->GetGeometryRef();

            if (geom != NULL && wkbFlatten(geom->getGeometryType()) == wkbPolygon) {

                OGRPolygon *OGRPoly = (OGRPolygon *) geom;

                auto start_extra_operation = clock::now();
                BoostPolygon poly = convertOGRPolygonToBoost(OGRPoly); // convert to BoostPolygon to use covered_by
                total_extra_opreations_time += (clock::now() - start_extra_operation);

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
                auto start_netcdf_query = clock::now();
                raster.getVar(start, count, subset.data());
                total_netcdf_queries_time += (clock::now() - start_netcdf_query);

                netcdfIOCount++;

                // Define default MAX value
                optional<double> currentMax = nullopt;
                // double currentMax = numeric_limits<double>::lowest();

                // Cell-by-cell Comparisons
                for (size_t i = 0; i < count[0]; i++) {
                    for (size_t j = 0; j < count[1]; j++) {

                        // Needed conversions
                        start_extra_operation = clock::now();
                        OGRPolygon cellPolygon = fromCellToPolygon(i + iymin, j + ixmin);
                        BoostPolygon boostedCellPolygon = convertOGRPolygonToBoost(&cellPolygon);
                        total_extra_opreations_time += (clock::now() - start_extra_operation);

                        if (bg::covered_by(boostedCellPolygon, poly)) {
                            double currentValue = subset[i * count[1] + j];
                            if (!currentMax.has_value() || currentValue > currentMax.value()) {
                                currentMax = currentValue;
                            }
                        }
                    }
                }

                start_extra_operation = clock::now();
                // Save pair <PolygonID, Max> only if a maximum value was found
                if (currentMax.has_value()) {
                    results[feature->GetFID()] = currentMax.value();
                } else {
                    // Exclude the current Polygon
                    numberOfExcludedPolygons++;
                    excludedPolygonsIDs.push_back(feature->GetFID());
                }
                total_extra_opreations_time += (clock::now() - start_extra_operation);
            }
            // Destroy feature
            OGRFeature::DestroyFeature(feature);
        }
        GDALClose(poDS);

        // Calculate execution time
        total_execution_time = clock::now() - start_total_time;

        // ================================     PRINT RESULTS    ================================
        
        printResultsTable(results);
        cout << "Number of polygons excluded: " << numberOfExcludedPolygons << endl;
        cout << "Polygons excluded: ";
        for (const long& polygonID : excludedPolygonsIDs) {
            cout << polygonID << " ";
        }

        cout << endl;
        cout << endl << "RESULTS OF TIME MEASURING" << endl << endl;
        cout << "Total execution time: " << total_execution_time.count() << " s "<< endl;
        cout << "Total NetCDF I/O operations: " << netcdfIOCount << endl;
        cout << "Total NetCDF I/O time: " << total_netcdf_queries_time.count() << " s " << endl;
        cout << "Total GDAL I/O operations: " << gdalIOCount << endl;
        cout << "Total GDAL I/O time: " << total_gdal_queries_time.count() << " s " << endl;
        cout << "Total extra operations time*: " << total_extra_opreations_time.count() << " s " << endl;
        cout << "------------------------------------------------------------------------" << endl;
        cout << "*Extra operations refers to other types of operations that are not related to NAIVE algorithm as such, for example: converting object types or collecting results." << endl;

        return 0;

    } catch (NcException& e) {
        cerr << "Error: " << e.what() << endl;
        return 1;
    }
}
