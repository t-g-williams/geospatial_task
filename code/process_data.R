### Geospatial task
### Import and do a basic descriptive analysis of climate, soil, and landcover data
### Author: Tim Williams

library(terra)
library(dplyr)
library(ggplot2)
library(rnaturalearth)
library(ncdf4)
library(sf)
library(rgdal)
library(data.table)
library(gridExtra)
library(cowplot)
library(ggjoy)
library(RColorBrewer)
library(pwr)

w_dir <- 'D:/My Drive/job_applications/2023_industry/climate_farmers/geospatial_task/'
setwd(w_dir)

main <- function() {
    ########## 1. Import and format the raw data ##########
    ## 1a. Select the region for the analysis --> North Island of New Zealand
    rgn <- select_NI_new_zealand()
    
    ## 1b. Biophysical data: load and crop (downloaded using python script "download_ERA5.py")
    era5_raw <- rast('data/biophysical/ERA5_download.nc')
    bio <- mask(crop(era5_raw, rgn), rgn)
    
    ## 1c. Landcover: crop the global data to the region
    lc_full <- rast('data/landcover/C3S-LC-L4-LCCS-Map-300m-P1Y-2020-v2.1.1.nc')
    lc <- mask(crop(lc_full, rgn), rgn)
    
    ## 1d. Soil: merge the files downloaded from SoilGrids
    soil <- merge_soilgrids_tifs(rgn)
    
    # write the data objects to file
    writeRaster(lc, 'data/landcover/landcover_region.tiff', overwrite=T)
    writeRaster(bio, 'data/biophysical/ERA5_download.tiff', overwrite=T)
    writeRaster(soil, 'data/soil/soil_combined.tif', overwrite=T)
    
    # make plots of the soil and landcover rasters
    plot_static_rasters(lc, soil)
    rm(lc_full, era5_raw) # clear some memory
    

    
    ########## 2. Align the data to a consistent spatial resolution ##########
    ### the biophysical data ("bio") has the lowest resolution
    ### --> resample the other data to this resolution
    
    ## 2a. landcover
    ## here we make an assumption! select the DOMINANT landcover cateogry in each aggregated pixel
    lc2 <- resample(lc, bio, method='mode') # "mode" : select the most common
    lc2 <- mask(lc2, bio[[1]]) # remove values that have no biophysical data
    
    ## 2b. soil
    ## as these data are continuous, we can take the mean within each pixel
    ## (note that this includes SOC values in non-agricultural land)
    soil2 <- resample(soil, bio, method='average')
    soil2 <- mask(soil2, bio[[1]]) # remove values that have no biophysical data
    
    # 2c. write these to file
    writeRaster(lc2, 'outputs/lc_resampled.tiff', overwrite=T)
    writeRaster(soil2, 'outputs/soil_resampled.tiff', overwrite=T)
    
    # 2d. finally, combine all into a single raster object
    ras <- c(lc2, soil2, bio) # using terra's stack function
    # add an ID variable layer
    tmp <- ras[[1]]
    names(tmp)[1] <- 'id'
    values(tmp) <- 1:size(tmp)
    values(tmp)[is.na(values(ras[['lccs_class']]))] <- NA # mask to the region
    ras <- c(ras, tmp)
    # write to file
    writeRaster(ras, 'outputs/combined_data.tiff', overwrite=T)
    
    
    
    ########## 3. Analysis and visualisation ##########
    # first, convert to long dataframe format (separately for static and temporal data)
    dfs <- raster_to_long_dfs(ras)
    dfs <- convert_units(dfs)
    
    # for the plotting, exclude landcover classes that only have a few pixels
    table(dfs[['static']]$lc_full)
    table(dfs[['static']]$lccs_class)
    lc_to_rm <- c(100, 150, 190, 200, 210, NA) # bare, urban, mosaic (tree/shrub), sparse vegetation, water
    
    # make timeseries plots
    plot_climate_timeseries(dfs, lc_to_rm)
    
    # histograms of SOC by land class
    plot_SOC(dfs[['static']], lc_to_rm)
    
    # map of rainfall
    plot_annual_rainfall(ras, dfs)
    

    
    ########## 4. Sample size calculation ##########
    # use the fine-resolution soil data
    N <- calculate_sample_size(soil, ci_level=0.95, perc_diff=10)
    print(paste('required sample size is', N, 'in each group'))
    
    
    
    ########## 5. extra task - climate risk ##########
    # classify pixels as high risk if they have low SOC and high precipitation variability
    # and low-risk is the opposite (high SOC, low precipitation variability)
    climate_risk_analysis(ras, dfs)
}


select_NI_new_zealand <- function() {
  # use the rnaturalearth package to get the country borders
  all_countries <- ne_countries(scale = "medium", returnclass = "sf")
  rgn_raw <- filter(all_countries, name=='New Zealand')
  
  # Need to manually extract the polygon representing the north island - it is the 11th polygon
  rgn <- st_cast(rgn_raw, "POLYGON")[11,]
  writeVector(vect(rgn), 'data/region.gpkg', overwrite=T) # write it to file to check
  
  # get the extent of rgn and write it to file
  rgn_extent <- ext(rgn)
  print(rgn_extent)
  ext_df <- data.frame('xmin'=rgn_extent$xmin, 'xmax'=rgn_extent$xmax, 
                       'ymin'=rgn_extent$ymin, 'ymax'=rgn_extent$ymax, row.names=c(1))
  write.csv(ext_df, 'data/region_extent.csv', row.names = F)
  
  return(rgn) 
}


merge_soilgrids_tifs <- function(rgn) {
  # the SoilGrids website only permits downloads in 2x2 degree formats
  # so we have multiple filenames, each containing a raster with a different extent
  # we want to load all of these rasters
  # and then collate them into a single object
  fns <- paste0('data/soil/out (', seq(0,15), ').tif')
  soil <- rast(fns[1])
  for (i in 2:16) {
    soil <- mosaic(soil, rast(fns[i]), fun=mean)
  }
  
  # clip to the region
  soil <- mask(soil, rgn)
  # correctly name the data
  names(soil) <- 'SOC_stock'
  
  return(soil)
}

plot_static_rasters <- function(lc, soil) {
    # plot the landcover and soil raster objects
    # first, convert to dataframes
    lc_df <- as.data.frame(lc, xy=T)
    soil_df <- as.data.frame(soil, xy=T)
    
    # plot the landcover
    # first, label the landcover classes
    landcover_lookup <- read.csv('data/landcover/key.csv')
    match_ixs <- match(lc_df$lccs_class, landcover_lookup$lccs_class)
    lc_df[,c('lc_full','lc_group')] <- landcover_lookup[match_ixs,c('name','group')]
    # change some types to "other"
    counts <- table(lc_df$lccs_class)
    excl <- names(counts)[counts < 6000]
    lc_df$lc_full[which(lc_df$lccs_class%in%excl)] <- "Other"
    lc_df$lc_full <- as.factor(lc_df$lc_full)
    p <- ggplot() +  
        geom_tile(data=lc_df, aes(x=x, y=y, fill=lc_full), alpha=1) + 
        # scale_fill_discrete(name='Landcover class') +
        scale_fill_brewer(palette='Paired',name='Landcover class') +
        coord_equal() +
        theme_map() +
        theme(legend.position="bottom") +
        theme(legend.key.width=unit(0.4, "cm")) +
      theme(plot.background = element_rect(fill = 'white', colour = 'white'))
    ggsave('outputs/plots/landcover_map.png', p, width=12, height=9)
    
    # plot the soil
    p2 <- ggplot() +  
        geom_tile(data=soil_df, aes(x=x, y=y, fill=SOC_stock), alpha=1) + 
        scale_fill_viridis_c(name='SOC stock (t/ha)') +
        coord_equal() +
        theme_map() +
        theme(legend.position="bottom") +
        theme(legend.key.width=unit(2, "cm"))
    ggsave('outputs/plots/soil_map.png', p2, width=8, height=7)
}


raster_to_long_dfs <- function(ras) {
  # create two dataframes from the raster
  # ONE for the temporal data (climate)
  # where each row represents a single combination of:
  # [pixel_id, variable, value, year]
  # and one for the static data (soil, landcover)
  # where each row represents [id, lcc_class, SOC_stock]
  
  df <- as.data.table(ras)
  row.names(df) <- 1:nrow(df)

  ## 1. temporal data
  # process the climate variables
  # which are labelled in {var}_{month} format
  bio_vars <- c('t2m','pev','stl1','e','tp')
  month_ixs <- seq(1,23*12)
  start_year <- 2000
  # loop over the months
  # (probably not the most computationally efficient, but it does the job)
  data_all <- list()
  for (m in month_ixs) {
    # select the data for this month
    cnames <- paste0(bio_vars, '_', m) # climate variables for this month
    cols_m <- c('id',cnames)
    data_m <- df[, ..cols_m]
    colnames(data_m) <- c('id', bio_vars)
    
    # melt into long format
    data_mm <- melt(data_m, id.vars='id')
    
    # add the time index
    data_mm$year <- start_year + (m-1)*1/12
    
    # add to the overall list of dataframes
    data_all[[m]] <- data_mm
  }
  # combine to single dataframe
  df_temporal <- do.call('rbind', data_all)
  fwrite(df_temporal, 'outputs/temporal_data_long.csv', row.names=F)
  
  ## 2. static data
  # simply select the columns
  static_vars <- c('id','lccs_class','SOC_stock')
  df_static <- df[, ..static_vars]
  # rename the landcover variables to something more intuitive
  # first, fix a few values
  df_static$lccs_class[df_static$lccs_class==81] <- 80
  df_static$lccs_class[df_static$lccs_class==12] <- 120
  landcover_lookup <- read.csv('data/landcover/key.csv')
  match_ixs <- match(df_static$lccs_class, landcover_lookup$lccs_class)
  df_static[,c('lc_full','lc_group')] <- landcover_lookup[match_ixs,c('name','group')]
  fwrite(df_static, 'outputs/static_data_long.csv', row.names=F)
  
  return(list('static'=df_static, 'temporal'=df_temporal))
}


convert_units <- function(dfs) {
  # convert the units as following:
  dfi <- dfs[['temporal']]
  
  # temperature and soil temperature: kelvin  --> degrees celsius
  ix1 <- dfi$variable %in% c('t2m','stl1')
  dfi$value[ix1] <- dfi$value[ix1] - 273
  
  # precipitation, PET, actual ET: metres/day --> millimetres/day
  ix2 <- dfi$variable %in% c('tp','e','pev')
  dfi$value[ix2] <- 1000 * dfi$value[ix2]
  
  # the evaporation ones store evaporation as NEGATIVE values
  # convert to positive
  ix3 <- dfi$variable %in% c('e','pev')
  dfi$value[ix3] <- -1 * dfi$value[ix3]
  
  dfs[['temporal']] <- dfi
  
  return(dfs)
}


plot_climate_timeseries <- function(dfs, lc_to_rm) {
  ### display timeseries for each climate variable
  ### grouped by landcover class
  
  # first, merge the dataframes to get the soil&LC data with the climate
  df <- dfs[['temporal']]
  vars_merge <- c('lccs_class','lc_full','SOC_stock')
  ix_match <- match(df$id, dfs[['static']]$id)
  df[, vars_merge] <- dfs[['static']][ix_match, ..vars_merge]
  
  # remove a few landcover classes that are not very common
  df <- df[-which(df$lccs_class %in% lc_to_rm),]
  
  # now group by the landcover class at each time point
  df_by_lc <- df %>% group_by(lc_full, year, variable) %>% summarise(
    value=mean(value), .groups='keep'
  )
  
  ## plot the entire timeseries
  single_temporal_plot(df_by_lc, '')
  ## and for a few years, so we can see it betteryear
  single_temporal_plot(df_by_lc %>% filter(floor(year)>=2020 & floor(year)<2022), '_2020-22')
  
  
  ## rolling average
  # the precipitation data are quite erratic: calculate a rolling average over the past 365 days
  df_by_lc2 <- df_by_lc %>% filter(variable=='tp')
  df_by_lc2$rolling_avg <- NA
  for (i in 1:nrow(df_by_lc2)) {
    print(i)
    if (df_by_lc2$year[i]>=2001) {
      df_tmp <- df_by_lc2 %>% filter(lc_full==df_by_lc2$lc_full[i])
      df_by_lc2$rolling_avg[i] <- mean(df_tmp$value[which(df_tmp$year>=df_by_lc2$year[i]-1 & df_tmp$year<=df_by_lc2$year[i])])
    }
  }
  # plot the rolling average
  p <- ggplot(df_by_lc2, aes(x=year, y=265*rolling_avg, color=lc_full)) +
    geom_line(linewidth=0.75) + labs(x='', y='1-year rolling sum of precipitation (mm/yr)') + theme_bw() + 
    theme(legend.position='bottom') + 
    guides(color=guide_legend(ncol=3))
  p
  ggsave('outputs/plots/rainfall_rolling_average.png', width=7, height=3.8)
}


single_temporal_plot <- function(df_in, name_ext) {
  # create a plot for each climate variable
  # using the input data in df_in
  
  bio_vars <- c('t2m','e','tp') # exclude potential ET (pev) and soil temp (stl1)
  bio_pretty <- c('Temperature (C)','Evaporation (mm/d)','Precipitation (mm/d)')
  
  plots <- c()
  for (i in 1:length(bio_vars)) {
    df_bio <- df_in %>% filter(variable==bio_vars[i])
    plots[[i]] <- ggplot(df_bio, aes(x=year, y=value, color=lc_full)) +
      geom_line() + labs(x='', y=bio_pretty[i]) + theme_bw() + 
      theme(legend.position='none')
    
    # remove the legend unless it is the last of the bio_vars
    if (i == length(bio_vars)) {
      plots[[i]] <- plots[[i]] + theme(legend.position='bottom') + 
        guides(color=guide_legend(ncol=3))
    }
  }
  
  # combine the plots into a single figure
  plot_all <- plot_grid(plotlist=plots, align = "v", nrow = 3, rel_heights = c(1,1,1.3))
  plot_all
  ggsave(paste0('outputs/plots/climate_timeseries', name_ext, '.png'), plot_all, width=7, height=10)
  
  return()
}


plot_SOC <- function(df_static, lc_to_rm) {
  # plot the distribution of soil quality data in each landcover class
  df_plt <- df_static %>% filter(!(lccs_class %in% lc_to_rm))

  # plot the distribution of SOC in each landcover class
  p1 <- ggplot(df_plt, aes(y=SOC_stock, x=lc_full)) +
    geom_boxplot(position='identity', alpha=0.5, fill='grey66') +
    labs(x='Landcover', y='SOC stock (t/ha)') + theme_classic() +
    theme(legend.position='bottom') + 
    guides(fill=guide_legend(ncol=3)) +
    lims(y=c(0,1.1*max(df_plt$SOC_stock))) +
    geom_hline(yintercept=0, color='black')
  p1
  ggsave('outputs/plots/SOC_by_landcover.png', p1, width=10, height=6)
}


calculate_sample_size <- function(soil, ci_level, perc_diff) {
  # extract the values
  soc <- values(soil)[!is.na(values(soil))]
  # remove the zero values --> these are for water etc.
  soc <- soc[-which(soc==0)]
  
  # calculate the moments of the distribution
  png(filename='outputs/plots/SOC_histogram.png', width=400, height=250)
  hist(soc)
  dev.off()
  mu <- mean(soc)
  sigma <- sd(soc)
  
  
  # calculate required sample size - alternative method
  # don't use this one because it is not about DIFFERENCES in SOC
  Z <- qnorm(1 - ((1 - ci_level) / 2)) # e.g. 1.96 for 95%
  abs_MOE <- mu*perc_diff/100 # e.g. 10% of the mean
  N <- (Z*sigma/abs_MOE)^2
  N <- ceiling(N) # round up
  
  
  # calculate required sample size
  # using a power analysis
  # assume power = 0.8 (this is standard)
  cohen_d <- (mu*perc_diff/100) / sigma # Effect size: difference between means (10% of the mean) divided by s.d.
  N_req <- ceiling(pwr.t.test(d=cohen_d, sig.level=1-ci_level, power=0.8)$n)
  
  # plot the sample size for a range of effect sizes
  mean_diffs <- seq(1, 20)
  N2 <- rep(NA, length(mean_diffs))
  for (i in 1:length(mean_diffs)) {
    cohen_d <- mean_diffs[i] / sigma
    N2[i] <- ceiling(pwr.t.test(d=cohen_d, sig.level=1-ci_level, power=0.8)$n)
  }
  df2 <- data.frame('group_difference'=mean_diffs, 'N'=N2)
  p <- ggplot(df2, aes(x=group_difference, y=N, label=N)) + 
    geom_vline(xintercept=mu/perc_diff, color='black') +
    geom_label() +
    labs(x='Expected difference in SOC due to management practice',
         y='Required sample size in each group') +
    theme_classic()
  p
  ggsave('outputs/plots/sample_sizes.png', width=8, height=4)
  
  return(N_req)
}


plot_annual_rainfall <- function(ras, dfs) {
  ## create a map of the average annual rainfall for each pixel
  
  # calculate annual precipitation for each pixel
  # note: this assumes each month has the same number of days, which is of course a little inaccurate
  df_tmp <- dfs[['temporal']] %>% filter(variable=='tp')
  df_tmp$year_round <- floor(df_tmp$year) 
  df_yr <- df_tmp %>% group_by(id, year_round) %>% summarise(
    tot_precip = mean(value) * 365, .groups='keep'
  )
  
  # calculate average precipitation for each pixel
  df_yr_avg <- df_yr %>% group_by(id) %>% summarise(mean_precip=mean(tot_precip))
  
  # get xy coords
  id_rast <- as.data.frame(ras[['id']], xy=T)
  id_rast$mean_precip <- df_yr_avg$mean_precip[match(id_rast$id, df_yr_avg$id)]
  
  # plot
  p <- ggplot() +  
    geom_tile(data=id_rast, aes(x=x, y=y, fill=mean_precip), alpha=1) + 
    # scale_fill_manual(values=c("high"='red',"moderate"='grey',"low"='blue'), name='Risk level') +
    coord_equal() +
    theme_map() +
    theme(legend.position="bottom") +
    theme(legend.key.width=unit(2, "cm"))
  p
  ggsave('outputs/plots/rainfall_map.png', p, width=9, height=6) 
}


climate_risk_analysis <- function(ras, dfs) {
  ## calculate a measure of climate risk for each pixel
  ## as the variance of annual precipitation
  ## then overlay this with the SOC data
  
  # calculate annual precipitation for each pixel
  # note: this assumes each month has the same number of days, which is of course a little inaccurate
  df_tmp <- dfs[['temporal']] %>% filter(variable=='tp')
  df_tmp$year_round <- floor(df_tmp$year) 
  df_yr <- df_tmp %>% group_by(id, year_round) %>% summarise(
    tot_precip = mean(value) * 365, .groups='keep'
  )
  # histogram (all pixels, all years)
  p <- ggplot(df_yr, aes(x=tot_precip)) + 
    geom_histogram(fill='grey66', color='black') +
    theme_classic() + labs(x='Annual precipitation (mm)')
  ggsave('outputs/plots/annual_precipitation.png', p, width=6, height=3)
  # joyplot showing change over time
  p <- ggplot(df_yr, aes(x=tot_precip, y=factor(year_round))) + 
    geom_joy() + theme_classic() + labs(x='Pixel-level total annual precipitation (mm)',y='')
  ggsave('outputs/plots/annual_precipitation_changes.png', p, width=5, height=6)
  
  
  ## (this is not part of the analysis)
  ## calculate mean annual temperature for each pixel
  df_tmp <- dfs[['temporal']] %>% filter(variable=='t2m')
  df_tmp$year_round <- floor(df_tmp$year) 
  df_yr_temp <- df_tmp %>% group_by(id, year_round) %>% summarise(
    mean_temp = mean(value), .groups='keep'
  )
  # basic plot
  p <- ggplot(df_yr_temp, aes(x=mean_temp)) + 
    geom_histogram(fill='grey66', color='black') +
    theme_classic() + labs(x='Average temperature (C)')
  ggsave('outputs/plots/average_temperature.png', p, width=6, height=3)
  # joyplot showing change over time
  p <- ggplot(df_yr_temp, aes(x=mean_temp, y=factor(year_round))) + 
    geom_joy() + 
    theme_classic() + labs(x='Pixel-level mean annual temperature',y='')
  ggsave('outputs/plots/average_temperature_changes.png', p, width=5, height=9)
  
  
  ## calculate the variance of annual precipitation for each gridcell
  df_sd <- df_yr %>% group_by(id) %>% summarise(
    precip_var = sd(tot_precip)
  )
  # histogram of this
  p <- ggplot(df_sd, aes(x=precip_var)) + 
      geom_histogram(fill='grey66', color='black') +
      theme_classic() + labs(x='Std. dev. of annual precipitation (mm)')
  ggsave('outputs/plots/annual_precipitation_variance.png', p, width=6, height=3)
  
  
  ## attribute the variance to the spatial object
  # first, format the variance in the correct order
  id_rast <- values(ras[['id']])
  rain_std_ras <- rep(NA, length(id_rast))
  rain_std_ras[match(df_sd$id, id_rast)] <- df_sd$precip_var
  # turn into raster object
  tmp_ras <- ras[['id']]
  values(tmp_ras) <- rain_std_ras
  names(tmp_ras)[1] <- 'precip_std'
  # append to the raster with SOC
  ras2 <- c(ras[['SOC_stock']], tmp_ras)
  
  
  ## assess associations between SOC and climate risk
  ## classify regions with (high climate risk, low SOC) as high risk
  ## and regions with (low climate risk, high SOC) as low risk
  # merge the two datasets together
  soc <- dfs[['static']]
  df_sd$SOC_stock <- soc$SOC_stock[match(df_sd$id, soc$id)]
  # calculate the tertiles for each variable to inform the risk levels
  soc_tertiles <- quantile(df_sd$SOC_stock, probs=c(0.33,0.66))
  rain_tertiles <- quantile(df_sd$precip_var, probs=c(0.33,0.66))
  # calculate risk levels
  df_sd$risk_level <- 1 # moderate
  df_sd$risk_level[(df_sd$SOC_stock<soc_tertiles[1]) & (df_sd$precip_var>rain_tertiles[2])] <- 2 # high
  df_sd$risk_level[(df_sd$SOC_stock>soc_tertiles[2]) & (df_sd$precip_var<rain_tertiles[1])] <- 0 # low
  df_sd$risk_level_text <- car::recode(df_sd$risk_level, '0="low";1="moderate";2="high"')
  
  # scatterplot (by pixel)
  # add dashed lines to show the tertiles
  # and text in the corners saying high-risk / low-risk
  p <- ggplot(df_sd, aes(x=SOC_stock, y=precip_var, color=risk_level_text)) + geom_point() + 
    theme_classic() +
    labs(x='SOC stock (t/ha)', y='Std. dev. of annual precipitation (mm)') +
    geom_vline(xintercept=soc_tertiles, linetype='dashed', linewidth=0.33) +
    geom_hline(yintercept=rain_tertiles, linetype='dashed', linewidth=0.33) +
    annotate('text', x=min(df_sd$SOC_stock), y=max(df_sd$precip_var), 
             label='High risk', hjust=0, vjust=1, color='red', size=10) + 
    annotate('text', x=max(df_sd$SOC_stock), y=min(df_sd$precip_var), 
             label='Low risk', hjust=1, vjust=0, color='blue', size=10) +
    scale_color_manual(values=c("high"='red',"moderate"='grey',"low"='blue'), name='Risk level') +
    theme(legend.position='none')
  p
  ggsave('outputs/plots/risk_scatterplot.png', p, width=8, height=6)
  
  ## add the risk classifications to the raster
  # re-order for raster object
  risk_ras <- rep(NA, length(id_rast))
  risk_ras[match(df_sd$id, id_rast)] <- df_sd$risk_level
  # turn into raster object
  tmp_ras <- ras[['id']]
  values(tmp_ras) <- risk_ras
  names(tmp_ras)[1] <- 'risk_level'
  # append to the raster with SOC
  ras2 <- c(ras2, tmp_ras)
  writeRaster(ras2, 'outputs/risk_analysis.tif', overwrite=T)
  

  # plot the raster of risk levels
  # color the pixels by risk level (red=high, blue=low, grey=moderate)
  ras_df <- as.data.frame(ras2[['risk_level']], xy=T)
  ras_df$risk_level_text <- car::recode(ras_df$risk_level, '0="low";1="moderate";2="high"')
  p <- ggplot() +  
    geom_tile(data=ras_df, aes(x=x, y=y, fill=risk_level_text), alpha=1) + 
    scale_fill_manual(values=c("high"='red',"moderate"='grey',"low"='blue'), name='Risk level') +
    coord_equal() +
    theme_map() +
    theme(legend.position="bottom") +
    theme(legend.key.width=unit(2, "cm"))
  ggsave('outputs/plots/risk_map.png', p, width=6, height=7)
}