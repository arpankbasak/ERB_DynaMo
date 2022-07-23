from czifile import CziFile
import pandas as pd
from tifffile import tifffile
import xml.etree.ElementTree as ET
import os
import sys
import xmltodict

image_loc=str(sys.argv[1])
f = list(os.listdir(str(image_loc)))
f = [i for i in f if not "._" in i]
output="./live_stock/"

for k in f:
  img_loc = str(image_loc+"/"+str(k))
  print("Converting "+str(img_loc)+" ....!")
  obj = CziFile(img_loc)
  z_layer = obj.shape[3]
  channel = obj.shape[1]
  mode = 0o755
  loc_img = img_loc.split("/")[-1]
  directory = loc_img.replace("./", "").replace(".czi", "")
  path_to_save = str(output+directory)
  try: 
      os.mkdir(path_to_save)
      print("Directory '% s' created" % directory)
  except OSError as error: 
      print(error)

  for i in range(0,channel):
    for j in range(0,z_layer):
      block = obj.asarray()
      img = block[0,i,0,j]
      tifffile.imwrite(path_to_save+"/"+directory+"_c"+str(i)+"_z"+str(j)+".tif", img)

  print("Storing metadata "+str(img_loc)+" ....!")
  
  # Storing the metadata
  meta = ET.fromstring(str(obj.metadata()).replace("\n", ""))

  temp = xmltodict.parse(str(obj.metadata()).replace("\n", ""))
  temp = temp['ImageDocument']['Metadata']

  temp_list = list(temp["Information"].items())
  df = pd.DataFrame(temp_list)
  df = pd.DataFrame(list(temp['Experiment']['ExperimentBlocks']['AcquisitionBlock']['AcquisitionModeSetup'].items())).T
  df = df.rename(columns = df.iloc[0]).loc[1:]

  df["Micron_X"] = df["DimensionX"].astype("float64")*df["ScalingX"].astype("float64")*1e6
  df["Micron_Y"] = df["DimensionY"].astype("float64")*df["ScalingY"].astype("float64")*1e6
  df["Micron_Z"] = df["DimensionZ"].astype("float64")*df["ScalingZ"].astype("float64")*1e6
  df["scale_value"] = df["ScalingX"].astype("float64")*1e6 * df["ScalingY"].astype("float64")*1e6*df["ScalingZ"].astype("float64")*1e6

  df.T.to_csv(path_to_save+"/"+directory+".txt", sep = "\t")

  print(" ...DONE!")