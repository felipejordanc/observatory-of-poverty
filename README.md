# observatory-of-poverty

## Required libraries
* Google Earth Engine API
* Geopandas
* Rasterio

## Usage

1. Download imagery
⋅⋅This section will download large tiles encompasing groups of villages to a Google Drive folder of your choice. Each tile is a mozaic of the available imagery within that year, selecting the best images in terms of their average level of NDVI (hence, images captures in full growth season are prioratized). By default, the Teillet C smooth correction method is applied to each image before the mozaic, but that can be turned off in step 1.4.
⋅⋅⋅⋅1. Create a local folder in the root of the repository. Git will ignore this folder.
⋅⋅⋅⋅2. Create file local/villages.csv. Must include three columns: village, with the name or ID of each village; lon, with latitude in decimal form (e.g 6.176874); and lon, with longitude in decimal form (e.g -14.114580).
⋅⋅⋅⋅3. Run d0_assignBatch.py script. See script for details.
⋅⋅⋅⋅4. Run d1_downloadBatches.py script, using the villages_batch.csv file created by the previos script. See script for details.
⋅⋅⋅⋅5. The previous script will download the imagery to a folder in your Google Drive. Download the folder and put the imagery in a folder in the local folder (e.g local/images).

2. Retrieve village image from download imagery
