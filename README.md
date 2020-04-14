# observatory-of-poverty

## Required libraries
* Google Earth Engine API
* Geopandas
* Rasterio

## Usage

1. **Download imagery**: This section will download large tiles encompasing groups of villages to a Google Drive folder of your choice. Each tile is a mozaic of the available imagery within that year, selecting the best images in terms of their average level of NDVI (hence, images captures in full growth season are prioratized). By default, the Teillet C smooth correction method is applied to each image before the mozaic, but that can be turned off in fourth step.
    1. Create a local folder in the root of the repository. Git will ignore this folder.
    2. Create file local/villages.csv. Must include three columns: village, with the name or ID of each village; lon, with latitude in decimal form (e.g 6.176874); and lon, with longitude in decimal form (e.g -14.114580).
    3. Run d0_assignBatch.py script. See script for details.
    4. Run d1_downloadBatches.py script, using the villages_batch.csv file created by the previos script. See script for details.
    5. The previous script will download the imagery to a folder in your Google Drive. Download the folder and put the imagery in a folder in the local folder (e.g local/images). Note that the script will also generate a file local/villages_loc.csv, which indicates what tile (or batch) contains the village.

2. **Retrieve villages' images from tiles**: The result of step 1 is a folder with images for large tiles that contains several villages. To create an image for a particular village, use the class villageImage defined in the village_image.py script. To initialize the class, you need to provide the name or ID of the village (same used in the village.csv file of step 1) and the longitude and latitude of the village. The instance can be used retrieve a subset of the original bands (which can optionally be pansharpened using the panchromatic band if available) as a numpy array, or save it to a tif image. See the end of the village_image.py script for a demonstration of the class usage.
