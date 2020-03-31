'''
This file retrieves satelitte imageries centered at the given coords of villages from the
GEE repository, created a quality mozaic by prioritizing imagery with higher average levels of NDVI,
and saves the quality mozaic in a Google Drive folder

Required argumentes:
    - folder: Name of folder on user's google drive where images will be saved
    - villages: Path to file with village coords relative to local folder located in root folder of this repo.
                File must include three columns: village, with the village name; lon, with the longitude of the village;
                and lat, with the latitude of the village
    - year: Year for which satellite images are collected to create quality mozaic
    - sensor: Sensor used (L5 for landsat 5, L7 for landsat 7)
    - bands: Output bands. Either all or rgb for all bands and 3 bands respectively
Flag arguments (optional):
    Pansharpenning:
    --panshapen: If sensor L7 is selected, uses panchromatic band to pansharpen the image using HSV method
    --no-pansharpen: Do not panshapen the image
    DEFAULT IS TO PANSHARPEN IMAGE
    --topocorrection: Make topographic correction using Teillet C method
    --no-topoocorrection: Do not make topographic correction
    DEFAULT IS NO TOPOGRAPHIC CORRECTION

The output image is a 244x244 tif image with 3 soil-reflectance bands (Near infrared, red, and green) in the [0,1] interval

Author:        Felipe Jordan
Date created:  12/16/2019
Last modified: 01/02/2020
'''

# Import the Earth Engine Python Package:
import sys
import os
import ee
import math
import pandas as pd
import argparse

# Initializa ee library
ee.Initialize()

# Add arguments
parser = argparse.ArgumentParser()
parser.add_argument("folder", help="Google Drive folder")
parser.add_argument("villages", help="Path of file with villages relative to local folder")
parser.add_argument("year", help="Year for which satellite images are collected")
parser.add_argument("sensor", help="Sensor. Can take values L5 and L7 for TM Landsat 5 and TME Landsat 7 respectively")
parser.add_argument("bands", help="Bands. Options are all or rgb.")
parser.set_defaults(bands='all')
parser.add_argument('--pansharpen'   , help="Use panchromatic band if available to pansharpen image when bands = rgb (default)", dest='pansharpen', action='store_true')
parser.add_argument('--no-pansharpen', help="Do not use panchromatic band if available to pansharpen image", dest='pansharpen', action='store_false')
parser.set_defaults(pansharpen=True)
parser.add_argument('--topocorrection'   , help="Make topographic correction (default)", dest='topocorrection', action='store_true')
parser.add_argument('--no-topocorrection', help='Do not make topographic correction', dest='topocorrection', action='store_false')
parser.set_defaults(topocorrection=True)
args = parser.parse_args()

########################## Helper functions and main class ---------
def keepClear(region):
  '''
  Initializes a function that is passed to the map method in GEE
  The child function takes an image as an argument, and mask all pixels that either:
    - The FMASK QA band assigns to it a medium or high confidence of cloud
    - The FMASK QA band assigns to it a medium or high confidence of cloud shadow
    - The saturation QA band determines that at least one band is saturated
  Finally, the function also adds to the image as a property the fraction of valid pixels within the region of interest,
  which is passed to the function that intializes the child function
  '''
  def keepClear_child(image):
    # Select FMASK QA band and select bits associated to clouds and clouds cover
    qa  = ee.Image(image).select('pixel_qa')
    qa_clouds            = extractQABits(qa,6,7)
    qa_cloudsShadows     = extractQABits(qa,3,3)
    # Select Saturation QA band and select bite associated to saturation
    qa2 = ee.Image(image).select('radsat_qa')
    qa_saturation        = extractQABits(qa2,1,1)
    # Create mask where valid pixels have low confidence of clound and clound shadow and are not saturated
    mask                 = qa_clouds.lte(1).And(qa_cloudsShadows.eq(0)).And(qa_saturation.eq(0))
    # Claculate fraction of valid pixels in the region and return image with QI property with fraction of valid pixels
    valid = mask.reduceRegion(ee.Reducer.sum(),region).get('pixel_qa')
    tot   = mask.reduceRegion(ee.Reducer.count(),region).get('pixel_qa')
    return image.updateMask(mask).copyProperties(image).set({'QI':ee.Number(valid).divide(tot)})
  return keepClear_child

def extractQABits(qaBand, bitStart, bitEnd):
  '''
  From a QA band, this function extract the information from bit bitStart to bit bitEnd and return in
  its decimal representation.
  '''
  numBits = bitEnd - bitStart + 1
  qaBits = qaBand.rightShift(bitStart).mod(2**numBits)
  return qaBits

def addNDVI(region):
  '''
  Initializes a function that is passed to the map method in GEE
  The child function takes an image as an argument, calculates the NDVI for each pixel, set negative values to zero, and saves in
  the image a property called NDVI that stores the average NDVI within a given region.
  The region is passed to the child_function from the main function.
  '''
  def addNDVI_child(image):
    ndvi = ee.Image(image).normalizedDifference(['B4', 'B3']).rename('ndvi')
    ndvi = ndvi.where(ndvi.lt(0),0)
    average_ndvi = ndvi.reduceRegion(ee.Reducer.mean(),region).get('ndvi')
    return image.copyProperties(image).set({'NDVI':ee.Number(average_ndvi)})
  return addNDVI_child

def panSharpen(region):
  '''
  Initializes a function that is passed to the map method in GEE

  The child function selects three bands, transforms them to the HSV space, replaces the value band (intensity)
  with the TOA panchromatic band (B8), and converts back to the RGB space. Before using the panchromatic band, it is
  reescaled to match the mean and standard deviation of the original image value band in the HSV space. Note that
  the range of wavelength captured by the panchromatic band in Landsat 7 covers the green, red, and near infra red
  bands. Hence, it make sence to use these bands (B2, B3, and B4) for the pansharpened image.

  The parent function is initialize with a region where the calculation of the average and standard deviation of the
  panchromatic band and the intensity of the color image are calculated.
  '''
  def panSharpen_child(image):
    #Convert the RGB bands to the HSV color space.
    hsv = ee.Image.rgb(image.select('B4').divide(ee.Image.constant(10000)),
                       image.select('B3').divide(ee.Image.constant(10000)),
                       image.select('B2').divide(ee.Image.constant(10000))).rgbToHsv()
    #Get panchromatic band from TOA collection, and then stretch it to match mean and std of value band in SR image
    pan           = image.select('B8')
    value         = hsv.select('value')
    pan_average   = ee.Image.constant(pan.reduceRegion(ee.Reducer.mean(),region,15).get('B8'))
    pan_std       = ee.Image.constant(pan.reduceRegion(ee.Reducer.stdDev(),region,15).get('B8'))
    val_average   = ee.Image.constant(value.reduceRegion(ee.Reducer.mean(),region,30).get('value'))
    val_std       = ee.Image.constant(value.reduceRegion(ee.Reducer.stdDev(),region,30).get('value'))
    pan_norm      = pan.subtract(pan_average).divide(pan_std)
    pan_stretched = pan_norm.multiply(val_std).add(val_average)
    #Swap in the panchromatic band and convert back to RGB.
    sharpened = ee.Image.cat([hsv.select('hue'), hsv.select('saturation'), pan_stretched]).hsvToRgb()
    return ee.Image(sharpened.copyProperties(image)).updateMask(image.mask().select('B2'))
  return panSharpen_child

def meters2degrees(lat):
  '''
  For a latitude (argument lat), returns a list with the equivalent in degrees of
  1 meter in the x direction (longitude) and 1 meter in the y directions (latitude).
  '''
  lon2met = 111111*math.cos(math.pi*lat/180)
  lat2met = 111111
  return [1/lon2met,1/lat2met]

def TerrainCorrection(scale,n_bands,smooth=5):
    '''
    Initializes a function that is passed to the map method in GEE
    The child function takes an image, applies the Teillet C smooth correction as describred in
    Riaño et al (2003) (https://ieeexplore.ieee.org/abstract/document/1206729), and returns the
    topographically corrected image.
    This function required the following parameters to initialize the child function:
        - Scale: scale of the final image
        - n_bands: Number of bands of the initial and final image
    '''

    def TerrainCorrection_child(img):
        degree2radian = 0.0174533
        #Region from footprint
        region = ee.Geometry.Polygon(ee.Geometry(img.get('system:footprint')).coordinates())
        #Extract solar zenith and calculate incidence angle (i)
        #Load USGS/SRTMGL1_003 DEM
        terrain = ee.call('Terrain', ee.Image('USGS/SRTMGL1_003')).clip(region)
        #Extract slope in radians for each pixel in the image
        p = terrain.select(['slope']).multiply(degree2radian).tan().divide(smooth).atan()
        #Extract solar zenith angle from the image
        z = ee.Image(ee.Number(img.get('SOLAR_ZENITH_ANGLE')).multiply(degree2radian))
        #Extract solar azimuth from the image
        az = ee.Image(ee.Number(img.get('SOLAR_AZIMUTH_ANGLE')).multiply(degree2radian))
        #Extract aspect in radians for each pixel in the image
        o = terrain.select(['aspect']).multiply(degree2radian)
        cosao = (az.subtract(o)).cos() #cos(ϕa−ϕo)

        #Calculate the cosine of the local solar incidence for every pixel in the image in radians (cosi=cosp*cosz+sinp*sinz*cos(ϕa−ϕo)
        cosi = img.expression('((cosp*cosz) + ((sinp*sinz)*(cosao)))',{'cosp': p.cos(),'cosz': z.cos(),'sinp': p.sin(),'sinz': z.sin(),'az' : az,'o' : o,'cosao': cosao})

        # Create the image to apply the linear regression.The first band is a constant, the second band the insidence angle cosi, and the next bands are the response variables (bands to which correction is being applied)
        # y = a + b*cosi
        reg_img = ee.Image.cat(ee.Image(1).rename('a'),cosi,img)
        #specify the linear regression reducer
        lr_reducer = ee.Reducer.linearRegression(**{'numX': 2,'numY': n_bands})
        #fit the model
        fit = reg_img.reduceRegion(**{'reducer': lr_reducer,'geometry': region,'scale': scale,'maxPixels': 1e10})

        # Calculate C corrector for each band: constant over slope
        coeff_array = ee.Array(fit.get('coefficients'))
        int = ee.Array(coeff_array.toList().get(0))
        slo = ee.Array(coeff_array.toList().get(1))
        C = int.divide(slo)
        Cimg = ee.Image.constant(C.toList())

        #Making the correction
        newimg = img.expression('((img * ((cosz) + C))/(cosi + C))',{'img': img,'cosz': z.cos(),'cosi': cosi,'C': Cimg})
        return newimg.copyProperties(img)
    return TerrainCorrection_child

# Main class that create image from sensor
class makeImage():
  def __init__(self,village,coords,year,sensor,bands,pansharpen=True,topocorrection=True):
    collections = {'L5':{'SR':['LANDSAT/LT05/C01/T1_SR','LANDSAT/LT05/C01/T2_SR'],'TOA':['LANDSAT/LT05/C01/T1_TOA','LANDSAT/LT05/C01/T2_TOA']},
                   'L7':{'SR':['LANDSAT/LE07/C01/T1_SR','LANDSAT/LE07/C01/T1_SR'],'TOA':['LANDSAT/LE07/C01/T1_TOA','LANDSAT/LE07/C01/T2_TOA']}}
    final_bands        = {'L5':['B1','B2','B3','B4','B5','B6','B7'],
                          'L7':['B1','B2','B3','B4','B5','B6','B7','B8']}
    pansharpenSensors = ['L7']
    self.village = village
    self.coords  = coords
    self.year    = year
    self.sensor  = sensor
    self.pansharpen     = pansharpen and (self.sensor in pansharpenSensors)
    self.topocorrection = topocorrection
    self.bands = bands
    if self.pansharpen:
      self.scale=15
    else:
      self.scale=30
    dx,dy       = meters2degrees(coords[1])
    self.region = ee.Geometry.Rectangle(self.coords[0]-122*self.scale*dx,
                                        self.coords[1]-122*self.scale*dy,
                                        self.coords[0]+122*self.scale*dx,
                                        self.coords[1]+122*self.scale*dy)
    try:
      self.collection = collections[sensor]
    except:
      print("Sensor must be L5 or L7")
    try:
      self.final_bands = final_bands[sensor]
    except:
      print("Sensor must be L5 or L7")

  def collect_images(self):
    # Select images and keep only those for which more than 95% of pixels are valid
    collection = ee.ImageCollection(self.collection['SR'][0])
    for c in self.collection['SR'][1:]:
      collection = collection.merge(ee.ImageCollection(c))
    collection = collection.filterDate('{0}-01-01'.format(self.year),
                                       '{0}-12-31'.format(self.year)).filterBounds(self.region)
    collection = collection.map(keepClear(self.region)).filter(ee.Filter.gte('QI',0.95))#.filter(ee.Filter.lte('SOLAR_ZENITH_ANGLE',60))

    # Add TOA if pansharpen is true
    def addTOA(img):
      pan = self.toa_collection().filter(ee.Filter.eq('LANDSAT_PRODUCT_ID',img.get('LANDSAT_ID'))).select('B8').first()
      img = img.addBands(pan.multiply(ee.Image.constant(10000)))
      return img

    if self.pansharpen:
        collection = collection.map(addTOA)

    return collection

  def length_collection(self):
    return self.collect_images().size().getInfo()

  def toa_collection(self):
    # return TOA collection
    toa_collection = ee.ImageCollection(self.collection['TOA'][0])
    for c in self.collection['TOA'][1:]:
      toa_collection = toa_collection.merge(ee.ImageCollection(c))
    return toa_collection

  def to_rgb(self):
    if self.pansharpen:
      # Add NDVI and sort according to average NDVI in image
      image     = self.collect_images().map(addNDVI(self.region)).sort('NDVI',False)
      if self.topocorrection:
          # Select NIR, RED and GREEN bands and apply terrain correction
          image     = image.select(['B8','B4','B3','B2']).map(TerrainCorrection(self.scale,4))
      # Pan sharpen image and select first not null for each pixel
      image     = image.map(panSharpen(self.region)).reduce(ee.Reducer.firstNonNull())
      # Select bands
      image     = image.select(['red_first','green_first','blue_first'])
    else:
      image     = self.collect_images().map(addNDVI(self.region)).sort('NDVI',False)
      if self.topocorrection:
          image     = image.select(['B4','B3','B2']).map(TerrainCorrection(self.scale,3))
      image     = image.reduce(ee.Reducer.firstNonNull())
      image     = ee.Image.rgb(image.select('B4_first'),
                               image.select('B3_first'),
                               image.select('B2_first')).divide(ee.Image.constant(10000)).updateMask(image.mask().select('B2_first'))
    return ee.Image(image).clip(self.region)

  def to_allBands(self):
    image     = self.collect_images().map(addNDVI(self.region)).sort('NDVI',False)
    if self.topocorrection:
        image     = image.select(self.final_bands).map(TerrainCorrection(self.scale,len(self.final_bands)))
    image     = image.reduce(ee.Reducer.firstNonNull()).divide(ee.Image.constant(10000))
    return ee.Image(image).clip(self.region)

  def to_drive(self):
    if self.bands=='all':
      img = self.to_allBands()
    elif self.bands=='rbg':
      img = self.to_rgb()
    #Save image to google drive
    export_image_toDrive = ee.batch.Export.image.toDrive(**{
                            'image': img,
                            'folder':'images',
                            'description': '{0}__{1}__{2}__{3}__{4}{5}'.format(self.village,self.year,self.sensor,self.bands,int(self.pansharpen),int(self.topocorrection)),
                            'dimensions': '244x244'})
    ee.batch.Task.start(export_image_toDrive)

################################### Main function
def main():
  #import villages files
  try:
    villages = pd.read_csv(os.path.join(os.pardir,'local',args.villages))
  except:
    print("Wrong file or file path to villages")
    sys.exit(1)

  #Create and save each image file to folder in Google Drive
  for i,row in villages.iterrows():
      try:
          name = row['village'].strip().lower().replace(' ','_')
          lon  = float(row['lon'])
          lat  = float(row['lat'])
      except:
          print('Villages file must include columns named "village", "lon", and "lat"')
          sys.exit(1)
      image = makeImage(name,[lon,lat],int(args.year),args.sensor,args.bands,args.pansharpen,args.topocorrection)
      if image.length_collection()>0:
          print("Saving image for village {0}".format(name))
          image.to_drive()
      else:
          print("Empty collection in village {0}. No image saved".format(name))


if __name__=="__main__":
  main()
