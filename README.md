# observatory-of-poverty

## Required libraries
* Google Earth Engine API
* Geopandas
* Rasterio

## Usage

1. **Download imagery**: This section will download large tiles encompasing a batch of villages to a Google Drive folder of your choice. Each tile is a mozaic of the available imagery within that year, selecting the best images in terms of their average level of NDVI (hence, images captured in full growth season are prioratized).
    1. Create a local folder in the root of the repository. Git will ignore this folder.
    2. Create file local/villages.csv. Must include three columns: village, with the name or ID of villages; lon, with latitude in decimal form (e.g 6.176874); and lon, with longitude in decimal form (e.g -14.114580).
    3. Run the d0_assignBatch.py script. See script for details.
    4. Run the d1_downloadBatches.py script, using as input the villages_batch.csv file created by the previos script. See script for details.
    5. The previous script will download the imagery to a folder in your Google Drive. Download the folder and put the imagery in the local folder (e.g local/images). Note that the script will also generate a file local/villages_loc.csv, which links each village to its batch number.

2. **Retrieve villages' images from tiles**: The result of step 1 is a folder with large images that contain several villages. To create an image for a particular village, use the class villageImage defined in the village_image.py script. To initialize the class, you need to provide the name or ID of the village (same used in the village.csv file of step 1) and the longitude and latitude of the village. By default, the class will load the local/villages_loc.csv file to locate the tile that contains the village. The instance can be used to retrieve a subset of the original bands (which can optionally be pansharpened using the panchromatic band if available) for image size as a numpy array, or save it as a tif file. See the end of the village_image.py script for a demonstration of the class usage.
