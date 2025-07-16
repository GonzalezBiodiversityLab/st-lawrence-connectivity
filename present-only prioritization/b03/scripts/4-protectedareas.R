#### a254: QC connectivity analysis
#### Script by Bronwyn Rayfield

#### 1. Filter Protected Areas

############################################################################################################################
# This code:                                                                                                               #
# - Removes aquatic protected areas from AP and RMN datasets                                                               #
# - Produces a vector layer of selected protected areas from the AP dataset only                                           #
# - Produces a vector layer of selected protected areas from the AP and RMN datasets                                       #
############################################################################################################################

# Load constants, functions, etc
source("./b03/scripts/0-0-constants.R")

#### Workspace ####

# Spatial data - Names
      # Study Area
studyArea_shp_Name <- file.path(b03ProcessedMapsDir, "b03-buffer-studyarea.shp")
studyArea_tif_Name <- file.path(b03ProcessedMapsDir, "b03-buffer-studyarea-30m.tif")
studyAreaFocal_shp_Name <- file.path(b03ProcessedMapsDir, "b03-studyarea.shp")
studyAreaFocal_tif_Name <- file.path(b03ProcessedMapsDir, "b03-studyarea-30m.tif")

      # Protected Areas
AP_Name <- file.path(b03RawMapsDir, "AP_REG_S_20210824.shp")
RMN_Name <- file.path(b03RawMapsDir, "RMN_20210608.shp")

      # Land cover
landCover_Name <- file.path(b03ProcessedMapsDir, "b03-landcover-30m.tif")

#### Set up GRASS Mapset #####
# If mapset has already been created
#initGRASS(gisBase=gisBase, gisDbase=gisDbase, location='BTSL_30m', mapset='ProtectedAreas', override=TRUE)

# Create location and PERMANENT mapset
initGRASS(gisBase=gisBase, gisDbase=gisDbase, location='BTSL_30m', mapset='PERMANENT', override=TRUE)

# Set projection info
execGRASS("g.proj", georef=studyArea_tif_Name, flags="c")

# Initialize new mapset inheriting projection info
execGRASS("g.mapset", mapset = "ProtectedAreas", flags="c")

# Import data
execGRASS("v.in.ogr", input=studyArea_shp_Name, output="rawData_vStudyArea", flags=c("overwrite", "o"))
execGRASS('r.in.gdal', input=studyArea_tif_Name, output='rawData_rStudyArea', flags=c("overwrite", "o"))
execGRASS("v.in.ogr", input=studyAreaFocal_shp_Name, output="rawData_vStudyAreaFocal", flags=c("overwrite", "o"))
execGRASS("r.in.gdal", input=studyAreaFocal_tif_Name, output="rawData_rStudyAreaFocal", flags=c("overwrite", "o"))
execGRASS("v.in.ogr", input=AP_Name, output="rawData_AP", snap=30, flags=c("overwrite", "o"))
execGRASS("v.in.ogr", input=RMN_Name, output="rawData_RMN", snap=30, flags=c("overwrite", "o"))
execGRASS('r.in.gdal', input=landCover_Name, output='landCover', flags=c("overwrite", "o"))

# Set region to b03 Lowlands and buffer
execGRASS('g.region', vector='rawData_AP')

#### Remove aquatic protected areas ####
# AP
execGRASS('v.extract', input='rawData_AP', where="NOT DESIGNOM LIKE 'Aire de concentration d''oiseaux aquatiques' AND NOT DESIGNOM LIKE 'Parc marin' AND NOT DESIGNOM LIKE 'Habitat du rat musqu%' AND NOT DESIGNOM LIKE '%(aire de nidification et bande de protection 0-200 m)'", output='AP_terrestrial', 'overwrite')

# RMN
      # Create a layer of water only
execGRASS('g.copy', raster=c('landCover', 'water'))
execGRASS('r.null', map='water', setnull='0-699')
execGRASS('r.null', map='water', setnull='701-1000')

      # Count number of water cells within each protected area
execGRASS('v.rast.stats', map='rawData_RMN', raster='water', column_prefix='water', method='sum')
execGRASS('v.db.addcolumn', map='rawData_RMN', columns='water_number integer')
execGRASS('v.db.update', map='rawData_RMN', column='water_number', query_column = 'water_sum / 700')

      # Count total number of cells within each protected area
execGRASS('v.rast.stats', map='rawData_RMN', raster='landCover', column_prefix='total', method='number')

      # Calculate percentage of protected area that falls in water
execGRASS('v.db.addcolumn', map='rawData_RMN', columns='water_percentage double')
execGRASS('v.db.update', map='rawData_RMN', column='water_percentage', query_column = '100 * water_number / total_number')
execGRASS('v.db.update', map='rawData_RMN', column='water_percentage', value='0', where='water_number IS NULL')
          
      # Remove RMN protected areas that fall in water
execGRASS('v.extract', input='rawData_RMN', where=paste('water_percentage <', waterThreshold), output='RMN_terrestrial', 'overwrite')

#### Merge polygons ####
# AP
      # Create a merge ID column with the same value for all polygons
execGRASS('v.db.addcolumn', map='AP_terrestrial', columns='mergeID integer')
execGRASS('v.db.update', map='AP_terrestrial', column='mergeID', value='1')

      # Dissolve on merge ID column
execGRASS('v.dissolve', input='AP_terrestrial', column='mergeID', output='AP_merged', 'overwrite')

      # Multipart to singlepart
execGRASS('v.category', input='AP_merged', output='AP_inter', option='del', cat=-1)
execGRASS('v.category', input='AP_inter', output='AP_merged', option='add', step=1, 'overwrite')

# RMN
      # Create a merge ID column with the same value for all polygons
execGRASS('v.db.addcolumn', map='RMN_terrestrial', columns='mergeID integer')
execGRASS('v.db.update', map='RMN_terrestrial', column='mergeID', value='1')

      # Dissolve on merge ID column
execGRASS('v.dissolve', input='RMN_terrestrial', column='mergeID', output='RMN_merged', 'overwrite')

      # Multipart to singlepart
execGRASS('v.category', input='RMN_merged', output='RMN_inter', option='del', cat=-1)
execGRASS('v.category', input='RMN_inter', output='RMN_merged', option='add', step=1, 'overwrite')

#### Produce Protected Area layer #1: AP areas only ####
# Calculate area of new polygons
execGRASS('v.db.addtable', map='AP_merged')
execGRASS('v.db.addcolumn', map='AP_merged', columns='area_ha double')
execGRASS('v.to.db', map='AP_merged', option='area', columns='area_ha', units='hectares', 'overwrite')

# Remove polygons that are not of interest
      # Set NULL areas in study area to zero
execGRASS('r.null', map='rawData_rStudyArea', null=0)

      # Polygons that are fully outside of the buffer
            # Calculate number of cells within the buffer, for each polygon
execGRASS('v.rast.stats', map='AP_merged', raster='rawData_rStudyArea', column_prefix='buffer', method='number', 'c')

            # Retain polygons that have at least one cell within the buffer
execGRASS('v.extract', input='AP_merged', where='buffer_number > 0', output='AP_buffer', 'overwrite')

      # Polygons that are smaller than the studyAreaThreshold
execGRASS('v.extract', input='AP_buffer', where=paste('area_ha >=', studyAreaThreshold), output='AP_inter', 'overwrite')

      # Polygons touching the buffer but outside the study area that are smaller than the bufferThreshold
            # Count, for each polygon, the number of cells within the study area
execGRASS('v.rast.stats', map='AP_inter', raster='rawData_rStudyArea', column_prefix='studyArea', method='sum')

            # Convert to area
execGRASS('v.db.addcolumn', map='AP_inter', columns='studyArea_ha double')
execGRASS('v.db.update', map='AP_inter', column='studyArea_ha', query_column='studyArea_sum * 900/10000')

            # Retain polygons that have more than the studyAreaThreshold within the study area OR their total area >= bufferThreshold
execGRASS('v.extract', input='AP_inter', where=paste0('(studyArea_ha >= ', studyAreaThreshold, ') OR (area_ha >= ', bufferThreshold, ')'), output='AP_final', 'overwrite')

      # Format columns
execGRASS('v.db.addcolumn', map='AP_final', columns='studyArea integer')
execGRASS('v.db.update', map='AP_final', column='studyArea', value='0', where=paste0('studyArea_ha < ', studyAreaThreshold))
execGRASS('v.db.update', map='AP_final', column='studyArea', value='1', where=paste0('studyArea_ha >= ', studyAreaThreshold))

      # Drop unnecessary columns
execGRASS('v.db.dropcolumn', map='AP_final', columns='buffer_number')
execGRASS('v.db.dropcolumn', map='AP_final', columns='studyArea_sum')
execGRASS('v.db.dropcolumn', map='AP_final', columns='studyArea_ha')

      # Remove intermediate products
execGRASS('g.remove', type='vector', name='AP_inter', 'f')

# Crop to region
      # Set region to that of buffer
execGRASS('g.region', raster='rawData_rStudyArea')

      # Create vector of region extent
execGRASS('v.in.region', output='region')

      # Crop to region
execGRASS('v.overlay', ainput='AP_final', binput='region', operator='and', output='AP', 'overwrite')

# Format output
execGRASS('v.reclass', input='AP', output='AP_inter', column='a_cat', 'overwrite')
execGRASS('v.db.addtable', map='AP_inter')
execGRASS('v.db.join', map='AP_inter', column='cat', other_table='AP', other_column="a_cat", subset_columns=c('a_area_ha', 'a_studyArea'))
execGRASS('v.db.renamecolumn', map='AP_inter', column=c('a_area_ha', 'area_ha'))
execGRASS('v.db.renamecolumn', map='AP_inter', column=c('a_studyArea', 'studyArea'))

# Remove intermediate products
execGRASS('g.rename', vector=c('AP_inter', 'AP_final'), 'overwrite')
execGRASS('g.remove', type='vector', name='AP', 'f')

# Save
      # Save csv of # of patches
nPatches <- v.get.att('AP_final', "@") %>%
  nrow(.) %>%
  as.data.frame()
colnames(nPatches) <- "NumberOfPatches"
write.csv(nPatches, file.path(b03ProcessedTabularDir, "num-patches-ap-buffer.csv"), row.names = F)

#      # Save shapefile
#execGRASS('v.out.ogr', input='AP_final', output=file.path(b03ProcessedMapsDir, 'protected-areas-ap.shp'), format='ESRI_Shapefile', flags=c('m', 'overwrite'))

#### Produce Protected Area layer #2: AP and RMN areas ####
# Set region to all of Quebec again
execGRASS('g.region', vector='rawData_AP')

# Create layer with both AP and RMN areas
execGRASS('v.overlay', ainput='AP_merged', binput='RMN_merged', output='all_terrestrial', operator='or', flags=c('overwrite'))


# Merge all overlapping and adjacent areas
      # Create a merge ID column with the same value for all polygons
execGRASS('v.db.addcolumn', map='all_terrestrial', columns='mergeID integer')
execGRASS('v.db.update', map='all_terrestrial', column='mergeID', value='1')

      # Dissolve on merge ID column
execGRASS('v.dissolve', input='all_terrestrial', column='mergeID', output='all_merged', 'overwrite')

      # Multipart to singlepart
execGRASS('v.category', input='all_merged', output='all_inter', option='del', cat=-1)
execGRASS('v.category', input='all_inter', output='all_merged', option='add', step=1, 'overwrite')

# Calculate area of new polygons
execGRASS('v.db.addtable', map='all_merged')
execGRASS('v.db.addcolumn', map='all_merged', columns='area_ha double')
execGRASS('v.to.db', map='all_merged', option='area', columns='area_ha', units='hectares', 'overwrite')

# Remove polygons that are not of interest
      # Polygons that are fully outside of the buffer
            # Calculate number of cells within the buffer, for each polygon
execGRASS('v.rast.stats', map='all_merged', raster='rawData_rStudyArea', column_prefix='buffer', method='number')

            # Retain polygons that have at least one cell within the buffer
execGRASS('v.extract', input='all_merged', where='buffer_number > 0', output='all_buffer', 'overwrite')

      # Polygons that are smaller than the studyAreaThreshold
execGRASS('v.extract', input='all_buffer', where=paste('area_ha >=', studyAreaThreshold), output='all_inter', 'overwrite')

      # Polygons touching the buffer but outside the study area that are smaller than the bufferThreshold
            # Count, for each polygon, the number of cells within the study area
execGRASS('v.rast.stats', map='all_inter', raster='rawData_rStudyArea', column_prefix='studyArea', method='sum')

            # Convert to area
execGRASS('v.db.addcolumn', map='all_inter', columns='studyArea_ha double')
execGRASS('v.db.update', map='all_inter', column='studyArea_ha', query_column='studyArea_sum * 900/10000')

            # Retain polygons that have more than the studyAreaThreshold within the study area OR their total area >= bufferThreshold
execGRASS('v.extract', input='all_inter', where=paste0('(studyArea_ha >= ', studyAreaThreshold, ') OR (area_ha >= ', bufferThreshold, ')'), output='all_final', 'overwrite')

      # Format columns
execGRASS('v.db.addcolumn', map='all_final', columns='studyArea integer')
execGRASS('v.db.update', map='all_final', column='studyArea', value='0', where=paste0('studyArea_ha < ', studyAreaThreshold))
execGRASS('v.db.update', map='all_final', column='studyArea', value='1', where=paste0('studyArea_ha >= ', studyAreaThreshold))

      # Drop unnecessary columns
execGRASS('v.db.dropcolumn', map='all_final', columns='buffer_number')
execGRASS('v.db.dropcolumn', map='all_final', columns='studyArea_sum')
execGRASS('v.db.dropcolumn', map='all_final', columns='studyArea_ha')

      # Remove intermediate products
execGRASS('g.remove', type='vector', name='all_inter', 'f')

# Crop to region
      # Set region to that of buffer
execGRASS('g.region', raster='rawData_rStudyArea')

      # Crop to region
execGRASS('v.overlay', ainput='all_final', binput='region', operator='and', output='all1', 'overwrite')

# Format output
execGRASS('v.reclass', input='all1', output='all_inter', column='a_cat', 'overwrite')
execGRASS('v.db.addtable', map='all_inter')
execGRASS('v.db.join', map='all_inter', column='cat', other_table='all1', other_column="a_cat", subset_columns=c('a_area_ha', 'a_studyArea'))
execGRASS('v.db.renamecolumn', map='all_inter', column=c('a_area_ha', 'area_ha'))
execGRASS('v.db.renamecolumn', map='all_inter', column=c('a_studyArea', 'studyArea'))

# Remove intermediate products
execGRASS('g.rename', vector=c('all_inter', 'all_final'), 'overwrite')
execGRASS('g.remove', type='vector', name='all1', 'f')

# Save
      # Save csv of # of patches
nPatches <- v.get.att('all_final', "@") %>%
  nrow(.) %>%
  as.data.frame()
colnames(nPatches) <- "NumberOfPatches"
write.csv(nPatches, file.path(b03ProcessedTabularDir, "num-patches-ap-rmn-multipart-buffer.csv"), row.names = F)

      # Save shapefile
#execGRASS('v.out.ogr', input='all_final', output=file.path(b03ProcessedMapsDir, 'protected-areas-ap-rmn-multipart.shp'), format='ESRI_Shapefile', flags=c('m', 'overwrite'))

# Multipart to singlepart
execGRASS('v.category', input='all_final', output='all_final_singlepart_inter', option='del', cat=-1)
execGRASS('v.category', input='all_final_singlepart_inter', output='all_final_singlepart', option='add', step=1, 'overwrite')
execGRASS('g.remove', type='vector', name='all_final_singlepart_inter', 'f')

# Save shapefile
execGRASS('v.out.ogr', input='all_final_singlepart', output=file.path(b03ProcessedMapsDir, 'protected-areas-150ha-btsl-900ha-buffer.shp'), format='ESRI_Shapefile', flags=c('m', 'overwrite'))


# Crop to within BTSL
execGRASS('v.overlay', ainput='all_final', binput='rawData_vStudyArea', operator='and', output='all_final_btsl', 'overwrite')

# Remove Gatineau Park
      # Isolate gatineau park (cat = 327)
execGRASS('v.extract', input='all_final', where='cat = 327', output='gatineauPark', 'overwrite')
      # Remove gatineau park from protected areas within the Lowlands
execGRASS('v.overlay', ainput='all_final_btsl', binput='gatineauPark', operator='not', output='all_final_btsl_no_gatineau', 'overwrite')

# Drop unnecessary columns
execGRASS('v.db.dropcolumn', map='all_final_btsl_no_gatineau', columns='a_cat')
execGRASS('v.db.dropcolumn', map='all_final_btsl_no_gatineau', columns='a_a_cat')
execGRASS('v.db.dropcolumn', map='all_final_btsl_no_gatineau', columns='a_a_area_ha')
execGRASS('v.db.dropcolumn', map='all_final_btsl_no_gatineau', columns='a_a_studyArea')
execGRASS('v.db.dropcolumn', map='all_final_btsl_no_gatineau', columns='a_b_cat')
execGRASS('v.db.dropcolumn', map='all_final_btsl_no_gatineau', columns='a_b_ID_NIV_01')
execGRASS('v.db.dropcolumn', map='all_final_btsl_no_gatineau', columns='a_b_NOM_PROV_N')
execGRASS('v.db.dropcolumn', map='all_final_btsl_no_gatineau', columns='a_b_SHAPE_Leng')
execGRASS('v.db.dropcolumn', map='all_final_btsl_no_gatineau', columns='a_b_SHAPE_Area')
execGRASS('v.db.dropcolumn', map='all_final_btsl_no_gatineau', columns='a_b_P_TERRE')
execGRASS('v.db.dropcolumn', map='all_final_btsl_no_gatineau', columns='a_b_sup_km2')
execGRASS('v.db.dropcolumn', map='all_final_btsl_no_gatineau', columns='b_cat')
execGRASS('v.db.dropcolumn', map='all_final_btsl_no_gatineau', columns='b_area_ha')
execGRASS('v.db.dropcolumn', map='all_final_btsl_no_gatineau', columns='b_studyArea')

# Save shapefile
execGRASS('v.out.ogr', input='all_final_btsl_no_gatineau', output=file.path(b03ProcessedMapsDir, 'protected-areas-150ha-btsl.shp'), format='ESRI_Shapefile', flags=c('m', 'overwrite'))


# Crop to within BTSL Focal (b03)
execGRASS('v.overlay', ainput='all_final_btsl_no_gatineau', binput='rawData_vStudyAreaFocal', operator='and', output='all_final_btsl_no_gatineau_focal', flags = c('t', 'overwrite'))
# Save shapefile
execGRASS('v.out.ogr', input='all_final_btsl_no_gatineau_focal', output=file.path(b03ProcessedMapsDir, 'protected-areas-150ha-btsl-focal.shp'), format='ESRI_Shapefile', flags=c('m', 'overwrite'))

# Rasterize PAs
# Set region to that of focal area
execGRASS('g.region', raster='rawData_rStudyAreaFocal')
# Rasterize
execGRASS("v.to.rast", input="all_final_btsl_no_gatineau_focal", use="val", value=1, output="all_final_btsl_no_gatineau_focal_30m", flags=c("overwrite"))
# Save raster
execGRASS('r.out.gdal', input='all_final_btsl_no_gatineau_focal_30m', output=file.path(b03ProcessedMapsDir, paste0('protected-areas-150ha-btsl-focal-', myResolution, 'm.tif')), format='GTiff', createopt='COMPRESS=LZW', flags=c('overwrite'))

