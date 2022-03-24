 

###################################################
#Libraries
#####################################################

list.of.packages <-
  c(
    "dplyr",
    "ggplot2",
    "readr",
    "lubridate",
    "data.table",
    "GGally",
    "here",
     "reshape2",
    "stringr",
    "conflicted",
    "lubridate",
    "naniar",
    "tibble",
    "tidyr"
  )

new.packages <-
  list.of.packages[!(list.of.packages %in% installed.packages()[, "Package"])]

#Download packages that are not already present in the library
if (length(new.packages))
  install.packages(new.packages)

if (length(new.packages))
  install.packages(new.packages, repos = c(CRAN="https://cran.r-project.org/"))

packages_load <-
  lapply(list.of.packages, require, character.only = TRUE)

#Print warning if there is a problem with installing/loading some of packages
if (any(as.numeric(packages_load) == 0)) {
  warning(paste("Package/s", paste(list.of.packages[packages_load != TRUE]), "not loaded!"))
} else {
  print("All packages were successfully loaded.")
}
conflict_prefer("here", "here")
conflict_prefer("map", "purrr")
conflict_prefer("select", "dplyr")
conflict_prefer("yday", "lubridate")
conflict_prefer("leap_year", "lubridate")
conflict_prefer("year", "lubridate")
conflict_prefer("filter", "dplyr")
conflict_prefer("yday", "lubridate")
conflict_prefer("month", "lubridate")

rm(packages_load, list.of.packages, new.packages)




#####################################################
# Find nearest stations
#####################################################

(df_loc <-
   read.csv(here::here("dat", "disease_loc.csv"), skip = 17) %>% 
   dplyr::select(-starts_with("X"))) 

(mtdt <- 
    read_csv( file = here("dat", "wth", "wth_loc_agwet.csv"))
)


# Computes distance using Haversine formula.
# Returns the result in meters.

haversine <- function( lat1, lon1, lat2, lon2, radius = 6371 ) {
  # Convert decimal degrees to radians
  lon1 = lon1 * pi / 180
  lon2 = lon2 * pi / 180
  lat1 = lat1 * pi / 180
  lat2 = lat2 * pi / 180
  
  # Haversine formula
  dlon = lon2 - lon1
  dlat = lat2 - lat1
  a = sin(dlat/2)**2 + cos(lat1) * cos(lat2) * sin(dlon/2)**2
  c = 2 * atan2(sqrt(a), sqrt(1-a))
  
  return( radius * c * 1000 )
}

#Add col with station name and col with distance to outbreaks
distances <- vector(mode = "numeric")

temp_st <- df_loc

for (i in seq_along(df_loc[["lat"]])) {
  
  # i = 2
  # y = 2
  for (y in seq_along(mtdt[["lat"]])) {
    distances[y]<-  round(haversine(df_loc[["lat"]][i], 
                                    df_loc[["lon"]][i],
                                    mtdt[["lat"]][y],
                                    mtdt[["lon"]][y] )/1000,2)
    
  }
  temp_st <- mtdt[order(distances),]
  
  df_loc$closest_1st[i] <- temp_st$stna[1]
  df_loc$dist_1st[i] <- distances[order(distances)][1]
  df_loc$closest_2nd[i] <- temp_st$stna[2]
  df_loc$dist_2nd[i]<- distances[order(distances)][2]
  df_loc$closest_3rd [i] <- temp_st$stna[3]
  df_loc$dist_3rd[i] <- distances[order(distances)][3]
  df_loc$closest_4th [i] <- temp_st$stna[4]
  df_loc$dist_4th[i] <- distances[order(distances)][4]
  df_loc$closest_5th [i] <- temp_st$stna[5]
  df_loc$dist_5th[i] <- distances[order(distances)][5]
  rm(i,y)
  
  
} 


df_loc$name <- as.character(df_loc$name)

write.csv(df_loc, here("dat", "dis_loc&station_disrances.csv"), row.names = F)

#####################################################
# Extract the data 
#####################################################

#Bio data: dis_long
load(file = here::here("dat", "raw", "final", "dis_long.RData"))

(dat_sum <-
    dis_long %>% 
    mutate(yr = year(date)) %>% 
    group_by(loc, yr) %>% 
    summarise() %>% 
    mutate(start = as.Date(paste(yr-1, "10", "01", sep = "-")),
           end = as.Date(paste(yr, "04", "30", sep = "-")))
)

(dis_long <- 
    dis_long %>% 
    select(-c(type, active)) %>% 
    mutate(yr = year(date)) %>% 
    unite(id, c("loc", "yr")))


# For each outbreak, subset the data from nearest station and add to 
# the new dataset
# Weather data
load(here("dat", "wth", "weather_hourly.RData"))
# wthh$yr <- year(wthh$datetime)
# wthh$mon <- lubridate::month(wthh$datetime)
wthh$date <- as.Date(wthh$datetime)

wthh %>% 
  group_by(stna) %>% 
  summarise(start = min(date),
            end = max(date))






wth_ls <- list()#list list for hourly wth
wthd_ls <- list()#list daily wth
 ch_ls<- list() #chill hour sums from different dates

   for (i in seq_along(dat_sum[["loc"]])) {
     # i=7
    #apothecia location and range of dates
    (loc <- dat_sum[i , "loc"] %>% pull())
    (dates <- seq.Date(pull( dat_sum[i , "start"]),
                       pull( dat_sum[i , "end"]), 
                       by = "day"))
    
    #vector of weather stations by distance
    (stnas <- 
        df_loc[df_loc$county == loc,grep("closest", names(df_loc))] %>% unlist())
    
    #Check if data from particular station is available within the given date range
    avaialble <- vector()
    for (y in seq(stnas)) {
      
      st_open <- wthh[ wthh$date %in% dates & wthh$stna == stnas[y] , "open"][1,]%>% pull()
      st_end <- wthh[ wthh$date %in% dates & wthh$stna == stnas[y] , "closed"][1,] %>% pull()
      
      wthh[  wthh$stna == stnas[y] , "closed"][1,]
      
      
      data_available <- 
        all(sapply(dates, function(x) x >= st_open && x <= st_end))
      avaialble [y] <- ifelse(data_available , 1, 0)
    }
    avaialble
    
    #select the closest station that has full data
    if(loc == "Island"){
      closest <- "Coupeville"
    }else{
      (closest <-  stnas[which(avaialble==1)[1]] %>% as.character())
    }
    
    
    wthf <- 
      wthh[wthh$stna == closest & wthh$date %in% dates,] 
    wthh[wthh$stna == closest ,] 
    
    
    
    #Distance of weather station used in the analysis
    # Find col with the closest station that has all data
    distance_col <- #some locations have the same name as the weather station
      if(length(which(df_loc[df_loc$county == loc, ] %in% closest))>1){
        which(df_loc[df_loc$county == loc, ] %in% closest)[2]
      }else{
        which(df_loc[df_loc$county == loc, ] %in% closest)
      }
    
    (wthf$dist <- df_loc[df_loc$county == loc,distance_col +1])
    
    
    #Add id for each loc/year
    wthf$id <- 
      dat_sum[i , c("loc", "yr") ] %>% paste(., collapse = "_")
    
    wthf <- 
      wthf %>% 
      mutate(doy = yday(date))
    
    #Mark which season 
    #different from year as data goes from October previous year
    wthf$season <- dat_sum[i , "yr"] %>% pull()
    
    
    #calculations for the initial observation and duration of presence of apothecia
    dis_sub <- 
      filter(dis_long, 
             id == unite(dat_sum[i , ], id, c("loc", "yr"))[ ,"id"] %>% pull())
    
    
    
    d_sum <- unite(dat_sum[i , ], id, c("loc", "yr"))[ ,"id"]
    
    d_sum$stna <-  closest
    d_sum$dist <- distance_col
    
    # Start of sporulation
    
    d_sum$cutoff <- "germ_start"
    wthf$germ_start <- 
      d_sum$date <-
      filter(dis_sub, stage == "germination") %>% 
      select(date) %>% 
      filter(row_number()==1 ) %>% pull()
    
    
    dff <- 
      filter(dis_sub,
             stage == "sporulation")
    
    if (any(dff$counts > 0)) {
      # Was there sporulation
      d_sum <- add_column(d_sum, sporulated = "yes", .after = "id")
      
      dates <-
        filter(dff, counts > 0) %>% select(date) %>% pull()
      
      # Sporulation
      spor_dur <-
        difftime(dates[length(dates)], dates[1]) %>%
        as.numeric(units = "days") + 1
      d_sum <- add_column(d_sum, spor_dur = spor_dur, .after = "sporulated")
      
      d_sum <-  mutate(d_sum[1,], 
                       cutoff ="spor_end",
                       date = dates[length(dates)]
      ) %>% bind_rows( d_sum, .) 
      
      
      
      d_sum <-  mutate(d_sum, 
                       cutoff ="spor_start",
                       date = dates[1]
                       
      ) %>% bind_rows(d_sum, .)
      
      d_sum <- 
        d_sum %>% 
        mutate( spor_start_prop = ifelse(leap_year(year(date)),
                                         yday(date) / (366 / 2),
                                         yday(date) / (365 / 2))
        )
      
      wthf$spore_start <- dates[1]
      wthf$spore_end <- dates[length(dates)]
      wthf$sporulated <- "yes"
      
    } else{
      # Was there sporulation
      d_sum <- add_column(d_sum, sporulated = "no", .after = "id")
      
      ( d_sum <- add_column(d_sum, spor_dur = 0, .after = "sporulated"))
      
      d_sum <-  mutate(d_sum[1,], 
                       cutoff ="spor_end",
                       date = as.Date(paste(dat_sum[i , "yr"] %>% pull(), "03", "15", sep = "-"))
      ) %>% bind_rows( d_sum, .) 
      
      d_sum <-  mutate(d_sum, 
                       cutoff ="spor_start",
                       date = as.Date(paste(dat_sum[i , "yr"] %>% pull(), "03", "15", sep = "-"))
                       
      ) %>% bind_rows(d_sum, .)
      
      
      
      
      d_sum <- 
        d_sum %>% 
        mutate( spor_start_prop = ifelse(leap_year(year(date)),
                                         yday(date) / (366 / 2),
                                         yday(date) / (365 / 2))
        ) 
      
      wthf$spore_start <- as.Date(paste(dat_sum[i , "yr"] %>% pull(), "03", "15", sep = "-"))
      wthf$spore_end <- as.Date(paste(dat_sum[i , "yr"] %>% pull(), "03", "15", sep = "-"))
      wthf$sporulated <- "no"
    }
    
    d_sum$cutoff <- factor(d_sum$cutoff, levels = c("germ_start","spor_start", "spor_end" ))
    
    d_sum <- distinct(d_sum)
    
    (  dis_sub <- 
        dis_sub %>% 
        group_by(stage) %>% 
        mutate(count_cum =cumsum(counts)))
    
    
    # dab - days after (chill hour) biofix
    lsdab <- split(wthf, wthf$date)
    for(z in seq(lsdab)) lsdab[[z]]$dab <- z
    wthf <- bind_rows(lsdab)
    
    
     
    
    ############################
    ### Units
    #############################
    
    wthf <-
      wthf %>%
      mutate(
        # Chilling units
        cu_ = ifelse(temp > 0 & temp < 7.2, 1, 0),
      ) 
    
    
    #Join apothecia development data
    wthf <-
      pivot_wider(dis_long, id_cols = c("id", "date"),
                  names_from =c( "stage"), 
                  values_from = "counts") %>% 
      filter(id ==  unique( wthf$id)) %>% 
      left_join(wthf, ., by = c("id", "date")) 
    
    
     # Define variable  containing all model names
    (levels_chill <- 
        select(wthf, 
               c(starts_with("cu_"),
               )) %>% colnames())
    
    
    # calculate daily values
    
    system.time(
      wthd <- 
        split(wthf, wthf$date, drop = TRUE) %>% 
        lapply(., function(x){
          # x <- split(wthf, wthf$date, drop = TRUE)[[1]]
          
          x
          
          y<- data.frame( date= date(x$datetime[1]),
                          stna =x$stna[1],
                          id =x$id[1],
                          sporulated = x$sporulated[1])
          y$temp <- mean(x$temp, na.rm = TRUE) %>% round(1)
          y$mintemp <- min(x$temp, na.rm = TRUE) %>% round(1)
          y$maxtemp <- max(x$temp, na.rm = TRUE) %>% round(1)
          
          y$soilt2in <- mean(x$soilt2in, na.rm = TRUE) %>% round(1)
          y$minsoilt2in <- min(x$soilt2in, na.rm = TRUE) %>% round(1)
          y$maxsoilt2in <- max(x$soilt2in, na.rm = TRUE) %>% round(1)
          
          y$soilm <- mean(x$soilm, na.rm = TRUE) %>% round(1)
          y$minsoilm <- min(x$soilm, na.rm = TRUE) %>% round(1)
          y$maxsoilm <- max(x$soilm, na.rm = TRUE) %>% round(1)
          
          y$soilt8in <- mean(x$soilt8in, na.rm = TRUE) %>% round(1)
          y$minsoilt8in <- min(x$soilt8in, na.rm = TRUE) %>% round(1)
          y$maxsoilt8in <- max(x$soilt8in, na.rm = TRUE, warning=FALSE) %>% round(1)
          
          y$rh <- mean(x$rh, na.rm = TRUE) %>% round(1)
          y$minrh <- min(x$rh, na.rm = TRUE) %>% round(1)
          y$maxrh <- max(x$rh, na.rm = TRUE, warning=FALSE) %>% round(1)
          
          y$wspd <- mean(x$wspd, na.rm = TRUE) %>% round(1)
          y$minwspd <- min(x$wspd, na.rm = TRUE) %>% round(1)
          y$maxwspd <- max(x$wspd, na.rm = TRUE, warning=FALSE) %>% round(1)
          
          y$lw <- sum(x$lw, na.rm = TRUE)
          y$solrad <- sum(x$solrad, na.rm = TRUE)
          y$rain <- sum(x$rain, na.rm = TRUE)
          
          y$dormant <- x$dormant[1]
          y$germination <- x$germination[1]
          y$emergence <- x$emergence[1]
          y$differentiation <- x$differentiation[1]
          y$sporulation <- x$sporulation[1]
          y$finished <- x$finished[1]
          
          
          y$frosthr <- 
            ifelse(any(x$temp<0), length(x$temp[x$temp<0]), 0)# hours temp less than zero
          
          
          # y$cu <-  sum(x$cu)
          # y$cu_zero <-  sum(x$cu_zero)
          
          y$germ_start <- x$germ_start[1]
          y$spore_start <- x$spore_start[1]
          y$spore_end <- x$spore_end[1]
          
           y <- 
            bind_cols(y,
                      select(x, all_of(levels_chill)) %>%
                        summarise_all(sum, na.rm = T))
          
          
          return(y)
        }) %>% bind_rows() 
    )
    
    #Clean up NaN and Inf 
    is.na(wthd) <- sapply(wthd, is.infinite)
    (wthd <- replace(wthd,is.na(wthd) , NA ))
    
    # rolling mean
    rlmean <- 10
    wthd[, paste("wtemp", rlmean, sep = "_")]  <-
      zoo::rollmean(wthd$temp, rlmean, mean, align = 'right', fill = NA)
    wthd[, paste("wmintemp", rlmean, sep = "_")]  <-
      zoo::rollmean(wthd$mintemp, rlmean, mean, align = 'left', fill = NA)
    wthd[, paste("wmaxtemp", rlmean, sep = "_")]  <-
      zoo::rollmean(wthd$maxtemp, rlmean, mean, align = 'left', fill = NA)
    
   
    
    
    # # Calculate Chill accumulations on different CHILL cutoff days
    chill_biofix <-
      c(
        paste("11", "01", sep = "-"),
        paste("11", "15", sep = "-"),
        paste("12", "01", sep = "-"),
        paste("12", "15", sep = "-"),
        paste("01", "01", sep = "-")
      )
    
    (  chill_dates <-
        c(as.Date(paste(wthf$season[1]-1, chill_biofix[1:4], sep = "-")),
          as.Date(paste(wthf$season[1], chill_biofix[5], sep = "-"))))
    
    chill_lss <- list()
    for (z in seq_along(chill_dates)) {
      chill <- chill_dates[[z]]
      chill_end <- d_sum[d_sum$cutoff == "spor_start", "date"] %>% unique() %>% pull()
      
      chill_lss[[z]] <-
        wthd %>%
        filter(date >= chill) %>%
        filter(date <= chill_end) %>%
        mutate(chill = chill) %>%
        mutate(across(starts_with("cu_"),  cumsum)) %>% 
        select(starts_with("cu_")) %>% 
        filter(row_number()==n()) %>% 
        bind_cols(d_sum[ d_sum$cutoff == "spor_start", ],.) %>%  
        add_column( chill_start = chill_biofix[z], .before = 1)
    }
    
    #chill units assessment 
    ch_ls[[i]] <-  
      bind_rows(chill_lss) %>% 
      mutate(chill_start = factor(chill_start, levels = chill_biofix))
    
    
    # Calculate the area under the curve for sporulation 
    aucdf <- dis_sub
    aucdf <- 
      aucdf[aucdf$stage == "sporulation", ] %>% 
      mutate(times = yday(date) - min(yday(date))) 
    ch_ls[[i]]["spor_int"] <-  sum(diff(aucdf$times)*zoo::rollmean(aucdf$prop,2))
    
    
    
    wth_ls[[i]] <- wthf
    wthd_ls[[i]] <- wthd
    dis_ls[[i]] <- dis_sub
    
    
      
    rm(wthf, dis_sub,lsdab, dff, loc, dates)
    print(i)
    
     
  }
 

 save(ch_ls, file = here("dat/ana/chill.RData"))
 save(wthd_ls, file = here("dat/ana/wthd_ls.RData"))

