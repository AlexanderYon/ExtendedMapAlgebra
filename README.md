-----

# Extended Map Algebra

## 🛠️ Dependencias del Sistema

Para compilar y ejecutar los programas, es necesario instalar las siguientes librerías de desarrollo:

  * `libgdal-dev`
  * `libnetcdf-c++4-dev`
  * `libnetcdf-dev`
  * `libboost-dev`

-----

## ⚙️ Programas Provistos

### 1\. `gen_raster.cpp` (Generador de Raster NetCDF)

Programa que genera un archivo raster en formato **NetCDF (`.nc`)** a partir de una especificación en texto plano.

#### Requisitos del Archivo de Entrada (`.txt`):

El formato del raster debe ser:

```
<n> <m>
a_11 a_12 ... a_1m
a_21 a_22 ... a_2m
 .    .   ...  .
 .    .   ...  .
a_n1 a_n2 ... a_nm
```

Donde:

  * `n` = número de filas.
  * `m` = número de columnas.
  * Los valores están separados por espacios y las filas por saltos de línea.

#### Compilación:

```bash
g++ gen_raster.cpp -o gen_raster -lnetcdf_c++4
```

#### Ejecución:

```bash
./gen_raster <archivo_raster_en_texto_plano> <archivo_raster_formato_netcdf>
```

-----

### 2\. `gen_polygons.cpp` (Generador de Shapefile)

Programa que genera un archivo **Shapefile** a partir de una especificación de polígonos en texto plano.

#### Requisitos del Archivo de Entrada (`.txt`):

El formato de los polígonos debe ser:

```
<n_poligonos>
<id0> <n_vertices>
<x0> <y0>
<x1> <y1>
...
<idn> <n_vertices>
...
```

#### Compilación:

```bash
g++ gen_polygons.cpp -I/usr/include/gdal -lgdal -o gen_polygons
```

#### Ejecución:

```bash
./gen_polygons <txtPath> <shpPath>
```

-----

### 3\. `naive.cpp` (Algoritmo Principal Ingenuo)

Programa principal que ejecuta el algoritmo ingenuo para encontrar el valor máximo del raster dentro de cada polígono.

#### Requisitos Previos:

Ejecutar `./gen_raster` y `./gen_polygons` **primero** para generar los archivos de datos (`.nc` y `.shp`).

#### Compilación:

```bash
g++ naive.cpp -o naive -lnetcdf_c++4 -lgdal -I/usr/include/gdal
```

#### Ejecución:

```bash
./naive <nc_raster> <shapefile>
```

-----

## ⚠️ Consideraciones Importantes

  * **Valores Máximos No Encontrados:** Si para algún polígono $P$, el programa `naive.cpp` no encuentra ninguna celda $C$ tal que $C$.`Within(P)` sea verdadero, el programa omitirá el polígono $P$ del reporte de resultados y lo agregará a la lista de polígonos excluídos.

-----

## 🔁 Ejemplo Flujo de Trabajo

Este es el orden recomendado para generar los datos y ejecutar el análisis:

```bash
# 1. Generar el archivo raster en formato NetCDF
./gen_raster raster_plain.txt raster.nc

# 2. Generar el archivo vectorial (Shapefile)
./gen_polygons polygons_plain.txt polygons

# 3. Ejecutar el algoritmo ingenuo
./naive raster.nc polygons/polygons.shp
```

La salida esperada del programa es la siguiente:

```bash
RESULTS
| Polygon ID          | Maximum              |
| id1                 | max1                 |
| id2                 | max2                 |
| id3                 | max3                 |
|  .                  |  .                   |
|  .                  |  .                   |
|  .                  |  .                   |
| idn                 | maxn                 |
Number of polygons excluded: <number>
Polygons excluded: <id0> <id1> <id2> ... <idk>

RESULTS OF TIME MEASURING

Total execution time: <time> s 
Total NetCDF I/O operations: <number>
Total NetCDF I/O time: <time> s 
Total GDAL I/O operations: <number>
Total GDAL I/O time: <time> s 
Total extra operations time*: <time> s 
------------------------------------------------------------------------

*Extra operations refers to other types of operations that are not related to NAIVE algorithm as such, 
for example: converting object types or collecting results.
```

Nótese que el tiempo de ejecución está medido en s