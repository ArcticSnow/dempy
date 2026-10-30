# GRASS

Install grass and launch grass console via `grass` command.

## Computing horizon for azimuth 45 [r.horizon](https://grass.osgeo.org/grass-stable/manuals/r.horizon.html)

```grass
r.import input=pooch/38e6e7f2f1f56999f3a217dcc7edff41-DEM_alp_L93_bilinear.tif output=dem
g.region raster=dem
r.horizon elevation=dem direction=45 step=0 output=horizon_45
r.out.gdal input=horizon_45 output=horizon45.tif format=GTiff
r.horizon elevation=dem direction=45 step=0 output=horizon_45
r.out.gdal input=horizon_45_045 output=horizon45.tif format=GTiff
```

## Computing sunmask


```grass
r.sunmask elevation=dem output=sunmask altitude=30 azimuth=165
```

- Calculation takes forever (still 0% after a few minutes)
- r.sunmask has been reported to be very slow by users, advice is to use r.sun instead.

## Sun irradiance

null values alc by sun irradiance correspond to shdowed areas.

```grass
r.sun elevation=dem day=355 time=10.5 incidout=incidence
r.mapcalc "ombres = if(isnull(incidence), 1, 0)"
r.out.gdal input=incidence output=incidence.tif format=GTiff
r.out.gdal input=ombres output=ombres.tif format=GTiff
```
