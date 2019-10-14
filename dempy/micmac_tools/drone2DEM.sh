# Script to compute DEM from Trond Eicken drone flights in Finse

module load micmac
module load gdal

# 1. Create a file "pic_name lat long elev"
mm3d GCPConvert AppInFile GCP_UTM_ellp.txt ChSys=SysUTM.xml@RTLFromExif.xml


#Use the GpsCoordinatesFromExif.txt file to create a xml orientation folder (Ori-RAWGNSS_N), and a file (FileImagesNeighbour.xml) detailing what image sees what other image (if camera is <50m away with option DN=50)
mm3d OriConvert "#F=N X Y Z" CAM_pos_UTM.txt RAWGNSS_N ChSys=SysUTM.xml@RTLFromExif.xml MTD1=1 NameCple=FileImagesNeighbour.xml DN=200   # in case coords in lat/long, DegreeWGS84@RTLFromExif.xml

#Find Tie points using 1/2 resolution image (best value for RGB bayer sensor)
mm3d Tapioca File FileImagesNeighbour.xml 2000

#filter TiePoints (better distribution, avoid clogging)
mm3d Schnaps .*JPG MoveBadImgs=1

mm3d Tapas FraserBasic .*JPG Out=Arbitrary SH=_mini

#Transform to  RTL system
mm3d CenterBascule .*JPG Arbitrary RAWGNSS_N Ground_Init_RTL

mm3d Campari .*JPG Ground_Init_RTL Ground_RTL EmGPS=[RAWGNSS_N,5] AllFree=1 SH=_mini

# mm3d SaisieAppuisPredicQT Ground_RTL GCP_UTM_ellp.xml GCP_picked.xml

# mm3d Campari .*JPG Ground_RTL Ground_GCP GCP=[GCP_picked.xml,0.3,GCPPos.xml,3] AllFree=1 SH=_mini


# mm3d ChgSysCo  .*JPG Ground_GCP RTLFromExif.xml@SysUTM.xml Ground_UTM

# mm3d OriExport Ori-Ground_UTM/O.*xml CameraPositionsUTM.txt AddF=1

# mm3d Malt Ortho .*JPG Ground_UTM ResolTerrain=0.1 ZoomF=2 EZA=1 DirMEC=MEC-GCP

# gdal_translate  -a_srs "+proj=utm +zone=32N +ellps=WGS84 +datum=WGS84 +units=m +no_defs" MEC-GCP/Z_Num9_DeZoom2_STD-MALT.tif OUTPUT/DEM_final.tif
# cp MEC-GCP/Z_Num9_DeZoom2_STD-MALT.tfw MEC-GCP/Correl_STD-MALT_Num_8.tfw
# gdal_translate  -a_srs "+proj=utm +zone=32N +ellps=WGS84 +datum=WGS84 +units=m +no_defs" MEC-GCP/Correl_STD-MALT_Num_8.tif OUTPUT/correl_final.tif

