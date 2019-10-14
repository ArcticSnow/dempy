# Document to apply artificial GCPs to unwarp drone DEM at Finse. 

# This pipeline is applied as a 2 processing steps. First compute an intial DEM, and then pick common points between a reference DEM and the initial DEM to finally compute an unwarped DEM 
# 	 - 0) copy the files in the project folder:
#			- sortGCPsIm.py
#			- sortGCPsIm.py
#			- unwarp_drone_part2.sh
#    - 1) an initial DEM must be computed using  DroneNadir.sh
#    - 2) pick common points (GCPs) in between a reference DEM and the first estimate DEM, store those in two shapefile with corresponding name. 
#    - 3) excecute this script.  ¤¤ adjust the input variables below ¤¤


# Prepare Wessel workspace
module load micmac
module load python/anaconda3
#module load gdal
source activate pybob

# IF NOT WORKING ON WESSEL:
#  - prepare a conda environment with python 3, gdal, geopandas, pandas, lxml, numpy
#  - activate this environment
#  - Have micmac installed

# Variable to check before launching
dem_ref=/home/arcticsnow/hycamp/Project/Finse/Finse-GIS/DEM/Dronefinse_DTM.tif		# path to the reference DEM (must be a GeoTif)
dem_target=OUTPUT/DEM_geotif.tif 		# path to the target DEM (must be a GeoTif)
GCP_ref_shapefile=GCP_ref.shp  			# Filename of the GCPs for the reference system (the final system)
GCP_target_shapefile=GCP_target.shp  	# Filename of the GCPs for the target system (the one you wanna change)

# Variables that can keep these default values as they are derived/local values 
GCP_ref_txt=GCP_ref.csv
GCP_target_txt=GCP_target.xyz
imPattern=*JPG
ori=Ori-Ground_UTM
GCP_name_prefix=p

# Add code to launch python script
python GCPs_sample_z.py -r $dem_target -s $GCP_target_shapefile -o $GCP_target_txt -gn $GCP_name_prefix
python GCPs_sample_z.py -r $dem_ref -s $GCP_ref_shapefile -o $GCP_ref_txt -gn $GCP_name_prefix

# 1. Convert GCPs xyz files to xml
mm3d GCPConvert AppInFile $GCP_ref_txt    # make sure the header #F=X Y N Z    or whatever the correct format is present at the top of the file

####################
# 2. Bob McNabb script adapted to our purpose. This section reproject the tiepoints found in each image into the new coordinate system.
cp -r $ori $ori-NoDist
tmp_autocal=$(ls $ori-NoDist/AutoCal*)
new_autocal=$(echo $tmp_autocal | sed 's_/_\\/_g')

sed -i 's/\(<CoeffDist.*>\)[^<>]*\(<\/CoeffDist.*\)/\10.0\2/' $tmp_autocal

for im in $imPattern; do
    echo "$im.xml"
    sed -i "s/\(<FileInterne.*>\)[^<>]*\(<\/FileInterne.*\)/\1$new_autocal\2/" $ori-NoDist/Orientation-$im.xml
    mm3d XYZ2Im $ori-NoDist/Orientation-$im.xml $GCP_target_txt NoDist-$im.txt
    mm3d XYZ2Im $ori/Orientation-$im.xml $GCP_target_txt Auto-$im.txt
done

# sortGCPsIm.py sorts through each images which GCPs are present withing the image or not, and spit out the xml file GCPPositionImage.xml. Python code initially written by Bob McNabb, adjusted to our case
python sortGCPsIm.py Auto-*.txt --no_distortion -o GCPPositionImage.xml -gn $GCP_name_prefix
rm -r $ori-NoDist
rm NoDist-*.txt
rm Auto-*.txt
####################

# 3. Compute transformation matrix in between reference and target coordinate system  WARNING: Change GCP_ref.xml based on GCP_ref_txt input!
mm3d GCPBascule .*JPG Ground_RTL Ground_GCP_Init GCP_ref.xml GCPPositionImage.xml

# 4. Re-estimate camera Orientation from this new transform (Ground_GCP_Init)
mm3d Campari .*JPG Ground_GCP_Init Ground_GCP GCP=[GCP_ref.xml,0.3,GCPPositionImage.xml,3] AllFree=1 SH=_mini

# 5. Compute the new DEM in the reference coordinate system
mm3d Malt Ortho .*JPG Ground_GCP ResolTerrain=0.1 ZoomF=2 EZA=1 DirMEC=MEC-GCP

# 6. Finish by running the gdal_translate command
gdal_translate  -a_srs "+proj=utm +zone=32N +ellps=WGS84 +datum=WGS84 +units=m +no_defs" MEC-GCP/Z_Num9_DeZoom2_STD-MALT.tif OUTPUT/DEM_final.tif
cp MEC-GCP/Z_Num9_DeZoom2_STD-MALT.tfw MEC-GCP/Correl_STD-MALT_Num_8.tfw
gdal_translate  -a_srs "+proj=utm +zone=32N +ellps=WGS84 +datum=WGS84 +units=m +no_defs" MEC-GCP/Correl_STD-MALT_Num_8.tif OUTPUT/correl_final.tif

# 7. Add code to delete tmp files and folders.