#####################################################
#Libraries
#####################################################

list.of.packages <-
  c(
    "here",
    "beepr",
    "ggplot2",
    "tibble",
    "tidyr",
    "data.table",
    "dplyr",
    "knitr",
    "imputeTS",
    "padr",
    "devtools",
    "readxl",
    "stringr",
    "lubridate",
    "readr",
    "zoo",
    "naniar",
    "conflicted",
    "weathermetrics"
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
conflict_prefer("year", "lubridate")
conflict_prefer("hour", "lubridate")
conflict_prefer("select", "dplyr")
conflict_prefer("month", "lubridate")
conflict_prefer("filter", "dplyr")
rm(packages_load, list.of.packages, new.packages)

################################################333
#Load and tidy data
#####################################################
my_files <- 
  list.files(here("dat", "wth", "raw2021"), pattern = "AWN_", recursive = TRUE)

my_files <-
  my_files[!my_files %in% stringr::str_subset(my_files, "Copy")]

my_files <- 
  paste0(here("dat", "wth", "raw2021"), "/", my_files)
(my_files <-  my_files[1:15])


ReadWeather <- function(path){
  read_csv(path, col_types = cols(
    TSTAMP_PST = col_character(),
    UNIT_ID = col_double(),
    STATION_NAME = col_character(),
    LATITUDE = col_double(),
    LONGITUDE = col_double(),
    ELEVATION_FEET = col_double(),
    AIR_TEMP_F = col_double(),
    SECOND_AIR_TEMP_F = col_double(),
    `RELATIVE_HUMIDITY_%` = col_double(),
    DEWPOINT_F = col_double(),
    LEAF_WETNESS = col_double(),
    PRECIP_INCHES = col_double(),
    WIND_DIRECTION_DEG = col_double(),
    WIND_SPEED_MPH = col_double(),
    SOLAR_RAD_WM2 = col_double(),
    SOIL_TEMP_2_IN_DEGREES_F = col_double(),
    SOIL_TEMP_8_IN_DEGREES_F = col_double(),
    `SOIL_MOIS_8_IN_%_VWC` = col_double(),
    AIR_PRESSURE_hPa = col_double()
  ))
}

all_csv <- lapply(my_files, function (x) {ReadWeather(x)})

wth <- 
 do.call("rbind", all_csv)                                   

rm(all_csv, my_files, cn, file_info)  


wth <-
  wth %>%
  rename(
    datetime = TSTAMP_PST,
    stno = UNIT_ID,
    stna = STATION_NAME ,
    lat = LATITUDE,
    lon = LONGITUDE,
    elev = ELEVATION_FEET, 
    temp_one = AIR_TEMP_F,
    temp_two = SECOND_AIR_TEMP_F,
    rh = `RELATIVE_HUMIDITY_%`,
    dwp = DEWPOINT_F,
    lw_resistance = LEAF_WETNESS,
    rain = PRECIP_INCHES,
    wdir = WIND_DIRECTION_DEG,
    wspd = WIND_SPEED_MPH,
    solrad = SOLAR_RAD_WM2,
    soilt2in = SOIL_TEMP_2_IN_DEGREES_F,
    soilt8in = SOIL_TEMP_8_IN_DEGREES_F,
    soilm = `SOIL_MOIS_8_IN_%_VWC`,
    airp = AIR_PRESSURE_hPa
  )

unique(wth$stna)

#Sort date and time columns
#Station MtVernon has a different time format
wth <- 
  split(wth, wth$stna) %>% 
  lapply(., function(x){
    if(x$stna[1]== "Mt.Vernon"){
      x$datetime <- mdy_hm(x$datetime)
    } else {
      x$datetime <- ymd_hms(x$datetime)
    }
    return(x)
  }) %>% bind_rows() %>% 
  dplyr::arrange(desc( stna))

summary(wth)

wth <- add_column(wth, yr = lubridate::year(wth$datetime), .after = "datetime")
wth <- add_column(wth, mon= lubridate::month(wth$datetime), .after = "yr")

#Replace 99999 with NA
wth[wth == 99999] <- NA

# drop all rows with NA values
wth <- filter(wth, rowSums(is.na(wth)) != ncol(wth))


#Extract metadata and save 
(mtdt <- 
wth %>% 
  group_by(stna) %>% 
  summarise(
    lat= unique(lat),
    lon = unique(lon),
    elev= unique(elev),
    open = min(as.Date(datetime)),
    closed = max(as.Date(datetime))) %>% 
  tbl_df()
  )
write_csv( mtdt,file = here("dat", "wth", "wth_loc_agwet.csv"))
rm(mtdt)

read_csv( file = here("dat", "wth", "wth_loc_agwet.csv"))

#reduce the size of data by removing the metadata
wth <-
wth%>%
  select(-c(stno,lat, lon, elev))

# Use the fist temp sensor. If there were no readings use the second
wth <- 
  wth%>%
  rowwise() %>% 
  mutate(temp = ifelse(!is.na(temp_one), temp_one,  temp_two)) %>% 
  select(-c(temp_one, temp_two)) %>% 
  ungroup()


#convert to metric system
wth$temp <- fahrenheit.to.celsius(wth$temp)
wth$soilt2in <- fahrenheit.to.celsius(wth$soilt2in)
wth$soilt8in <- fahrenheit.to.celsius(wth$soilt8in)
wth$rain <- 
  convert_precip(wth$rain,old_metric = "inches", new_metric = "mm",round = 2)

#convert solar radiation to MJ/m/



#Impute the missing data 
infil_gap <- 12 #Maximum length of the infill gap
wth$temp <-
  round(na.spline(wth$temp, na.rm = FALSE, maxgap = infil_gap), 1)

wth$dwp <-
  round(na.spline(wth$dwp, na.rm = FALSE, maxgap = infil_gap), 1)

wth$lw_resistance <- 
  na.approx(wth$lw_resistance,maxgap = 3)

wth$rain <- 
  na.approx(wth$rain,maxgap = 3)

wth$wdir <- 
  na.approx(wth$wdir,maxgap = 3)
wth$wspd <- 
  na.approx(wth$wspd,maxgap = 3)

wth$soilt2in <-
  round(na.spline(wth$soilt2in, na.rm = FALSE, maxgap = infil_gap), 1)
wth$soilt8in <-
  round(na.spline(wth$soilt8in, na.rm = FALSE, maxgap = infil_gap), 1)

wth$rh <-
  round(na.spline(wth$rh, na.rm = FALSE, maxgap = infil_gap), 0)
wth$rh  <- sapply(wth$rh, function(x)
  ifelse(x > 100, x <- 100, x))




save(wth,file =  here::here("dat", "wth", "weather.RData"))

################################
# Calculate the hourly averages 
###################################
load(here("dat", "wth", "weather.Rdata"))

# The period was treated as leaf was wet if the electrical resistance was >= 0.4. 
# These were then averaged over hourly periods. 
wth$lw <- ifelse(wth$lw_resistance > .35, 1,0)


wth %>% 
  summary()

#Soil temperatures below -50C are probably errors 
wth %>% 
  filter(soilt2in < (-50))


# temperatures below -50C are probably errors 
hist(wth$temp)
wth %>%  filter(temp < (-40))

wth %>% arrange(temp) %>% select(datetime, stna, temp) %>% head(1500) 



wth <- 
wth %>% 
  mutate(
    soilt2in = ifelse(soilt2in <(-50), NA, soilt2in),
    temp = ifelse(temp <(-50), NA, temp)
    )


# Hourly averages are calculated by reducing the time stamp for 1 second and 
# rounding up to the nearest higher unit
wth$datetime <- ceiling_date(wth$datetime-100, 
                     unit = "hour",change_on_boundary = TRUE)

wthh <- 
  wth %>% 
  group_by(stna, datetime) %>% 
  summarise(
    # yr = unique(yr),
    temp = mean(temp, na.rm = TRUE ),
    rh =  mean(rh, na.rm = TRUE),
    dwp =  mean(dwp, na.rm = TRUE),
    lw =  mean(lw, na.rm = TRUE),
    rain = sum(rain, na.rm = TRUE),
    wspd = mean(wspd, na.rm = TRUE),
    wdir = mean(wdir, na.rm = TRUE),
    solrad = sum(solrad, na.rm = TRUE),
    soilt2in = mean(soilt2in),
    soilt8in = mean(soilt8in),
    soilm = mean(soilm),
    airp = mean(airp, na.rm = TRUE)
  ) %>% 
  ungroup()


#Impute the missing data 
infil_gap <- 12 #Maximum length of the infill gap
wthh$temp <-
  round(na.spline(wthh$temp, na.rm = FALSE, maxgap = infil_gap), 1)

wthh$dwp <-
  round(na.spline(wthh$dwp, na.rm = FALSE, maxgap = infil_gap), 1)

wthh$lw <- 
  na.approx(wthh$lw,maxgap = 3)

wthh$rain <- 
  na.approx(wthh$rain,maxgap = 3)

wthh$wdir <- 
  na.approx(wthh$wdir,maxgap = 3)
wthh$wspd <- 
  na.approx(wthh$wspd,maxgap = 3)

wthh$soilt2in <-
  round(na.spline(wthh$soilt2in, na.rm = FALSE, maxgap = infil_gap), 1)
wthh$soilt8in <-
  round(na.spline(wthh$soilt8in, na.rm = FALSE, maxgap = infil_gap), 1)

wthh$rh <-
  round(na.spline(wthh$rh, na.rm = FALSE, maxgap = infil_gap), 0)
wthh$rh  <- sapply(wthh$rh, function(x)
  ifelse(x > 100, x <- 100, x))

wthh$rh  <- sapply(wthh$rh, function(x)
  ifelse(x < 0, NA, x))


wthh %>% 
  mutate(yr = lubridate::year(datetime)) %>% 
  dplyr::filter(yr>2010) %>%
  mutate(mon = lubridate::month(datetime)) %>% 
  dplyr::filter(mon %in% c(9:12,1:5)) %>% 
  select(., -c(airp,wspd,yr, datetime, dwp, mon)) %>% 
  gg_miss_fct(x = ., fct = stna) 

  ggsave(here::here("out", "wth", "na_post2010"), device = "png")

wthh %>% 
  mutate(yr = lubridate::year(datetime)) %>% 
  dplyr::filter(yr>2010) %>%
  mutate(mon = lubridate::month(datetime)) %>% 
  dplyr::filter(mon %in% c(9:12,1:5)) %>% 
  select(., -c(airp,wspd,yr, soilt2in, datetime, dwp, mon, soilm, wdir)) %>% 
  gg_miss_fct(x = ., fct = stna)


mtdt <- 
  read_csv( file = here("dat", "wth", "wth_loc_agwet.csv"))

wthh <-
  left_join(wthh, mtdt, by = "stna")

##############################
##### Leaf wetness checks

wthh %>% 
  sample_frac(.05) %>% 
  mutate(lw = ifelse(lw > .4, 1, 0)) %>% 
  ggplot(aes(rh,lw))+
  geom_point()+
  stat_smooth(method="glm",family="poisson",link="probit", formula=y~x,col="red")

dtt <- 
wthh %>% 
  sample_frac(.1) %>% 
  mutate(lw = ifelse(lw > .5, 1, 0)) 
  
my_fit <- 
  glm(lw ~ rh, data = dtt, na.action = na.exclude,
              family = "binomial")
summary(my_fit)
pred <- predict(my_fit, type = "response")
pred_df <- data.frame(rh = dtt$rh, lw =pred)
ggplot(dtt, aes(x = rh, y = lw)) +
  geom_point() +
  geom_smooth(method = "glm", 
              method.args = list(family = "binomial")) +
  geom_point(data = pred_df, aes(rh, lw), colour = "blue") +
  geom_hline(data = data.frame(c(0.25, 0.50, 0.75)),
             aes(yintercept = c(0.25, 0.50, 0.75)),
             colour = "darkgrey", linetype = "dashed")

my_fit <- 
  glm(lw ~ rain, data = dtt, na.action = na.exclude,
      family = "binomial")
summary(my_fit)
pred <- predict(my_fit, type = "response")
pred_df <- data.frame(rain = dtt$rain, lw =pred)
ggplot(dtt, aes(x = rain, y = lw)) +
  geom_point() +
  geom_smooth(method = "glm", 
              method.args = list(family = "binomial")) +
  geom_point(data = pred_df, aes(rain, lw), colour = "blue") +
  geom_hline(data = data.frame(c(0.25, 0.50, 0.75)),
             aes(yintercept = c(0.25, 0.50, 0.75)),
             colour = "darkgrey", linetype = "dashed")
dtt %>% 
  filter(rain<5) %>% 
  ggplot()+
  geom_histogram(aes(rain), bins = 100)+
  facet_wrap(~lw)



wthh %>% 
  mutate(yr = year(datetime)) %>% 
  group_by(stna) %>% 
  summarise(
    max_lw = max(lw, na.rm = T),
    mean_lw = mean(lw, na.rm = T),
  )
nrow(wthh[wthh$lw >1,] %>% distinct()) /nrow(wthh)




# Check for >18C ----------------------------------------------------------
wthh$wet_dur <- ifelse(wthh$rh >= 90 | wthh$rain > .1, 1, 0)



wthh %>%
  mutate(mn = lubridate::month(datetime)) %>%
  filter(mn %in% c(2:4))

wthh %>%
  mutate(mn = lubridate::month(datetime)) %>%
  filter(mn %in% c(2:4)) %>%
  filter(temp > 18 & wet_dur == 1)
 
summary(wthh)

save(wthh, file =  here::here("dat", "wth", "weather_hourly.RData"))


####################################
#Upscaling to daily resolution
####################################

load(here("dat", "wth", "weather_hourly.RData"))

wthh <- 
  wthh %>% 
  mutate(date = date(datetime)) %>% 
  tidyr::unite("id", c(date,stna),remove = FALSE )
  
  
wthls <- 
  split(wthh, wthh$id, drop = TRUE)

lenghts <- sapply(wthls, nrow) %>% unlist()
lenghts[lenghts<24] %>% length()
ss <- lenghts>3
wthls <- wthls[ss]

system.time(
  wthd <- 
lapply(wthls, function(x){
  # x <- wthls[[2333]]
  # x <- wthls[[33333]]
  # x <- wthls[[3]]
  # x <- wthls[[4]]
  x
  
  y<- data.frame( date= date(x$datetime[1]),
                  stna =x$stna[1])
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
  return(y)
}) %>% bind_rows() 
  )


#Clean up NaN and Inf 
is.na(wthd)<-sapply(wthd, is.infinite)
wthd <- replace(wthd,is.na(wthd) , NA )

# save(wthd,file =  here::here("dat", "wth", "weather_day.RData"))
load(file =  here::here("dat", "wth", "weather_day.RData"))


summary(wthd)
nrow(wthd)
nas <- 
wthd %>% 
  mutate(yr = year(date),
         mon = month(date)) %>% 
  dplyr::filter(yr>2014) %>% 
  group_by(stna) %>% 
  naniar::miss_summary() %>% 
  select(miss_case_table) %>% 
  tbl_df()
nas[[1]]

# full year
ttl <- "2015-2020_daily" 
wthd %>% 
  mutate(yr = year(date),
         mon = month(date)) %>% 
  dplyr::filter(yr>2014) %>%
  # dplyr::filter(mon %in% c(9:12, 1:5)) %>%
  select(., c(lw, rh, temp, soilt8in, soilt2in, solrad, rain, soilm, stna)) %>% 
  gg_miss_fct(x = ., fct = stna)+
  ggtitle(ttl)
ggsave(here::here("out","na", paste0("na_", ttl, ".png")))

ttl <- "2015-2020_low_na_daily" 
wthd %>% 
  mutate(yr = year(date),
         mon = month(date)) %>% 
  dplyr::filter(yr>2014) %>%
  # dplyr::filter(mon %in% c(9:12, 1:5)) %>%
  select(., c(lw, rh, temp, soilt8in, rain, stna)) %>% 
  gg_miss_fct(x = ., fct = stna)+
  ggtitle(ttl)
ggsave(here::here("out","na", paste0("na_", ttl, ".png")))

# From Oct to May
ttl <- "2015-2021_daily_season" 
wthh %>% 
  mutate(yr = year(datetime),
         mon = month(datetime)) %>% 
  dplyr::filter(yr>2014) %>%
  dplyr::filter(mon %in% c(9:12, 1:5)) %>%
  select(., c(rh, temp,rain,  soilt8in, soilt2in,  soilm,solrad,lw, stna)) %>%
  rename("Temperature (°C)" = temp,
         "Solar Radiation (MJ/cm2)" = solrad,
         "Soil Temperature (8'')" = soilt8in, 
         "Soil Temperature (2'')" =soilt2in,  
         "Soil Moisture (2'')" =soilm,
         # "Solar Radiation" = solrad,
         "Relative Humidity (%)" = rh,
         "Precipitation (l/m2)" = rain,
         
         "Leaf wetness" =lw
  ) %>%
  gg_miss_fct(x = ., fct = stna)+
  
  ylab("Weather Variable")+
  xlab("Weather Station")+
  theme(legend.position = "top")
ggsave(here::here("out","na", paste0("na_", ttl, ".png")),
       width = 8,
       height =4.5)

ttl <- "2015-2021_daily_season_low_na" 
wthd %>% 
  mutate(yr = year(date),
         mon = month(date)) %>% 
  dplyr::filter(yr>2014) %>%
  dplyr::filter(mon %in% c(9:12, 1:5)) %>%
  select(., c(lw, rh, temp, soilt8in, rain, stna)) %>%
  gg_miss_fct(x = ., fct = stna)+
  ggtitle(ttl)
ggsave(here::here("out","na",  paste0("na_", ttl, ".png")))


#Calculate summaries per variable, year and station
# browseURL("https://cran.r-project.org/web/packages/naniar/vignettes/naniar-visualisation.html")
library("naniar")

gg_miss_upset(wthh)

gg_miss_var(wthh,facet = stna)
gg_miss_var(select(wthh, -c(airp,wspd)),facet = stna)

gg_miss_fct(x = wthh, fct = stna)
gg_miss_fct(x = select(wthh, -c(airp,wspd)), fct = stna)


wthh %>%
  mutate(yr = year(date),
         mon = month(date)) %>% 
  
  group_by(stna, yr) %>%
  miss_var_summary()



png(here::here("out","wth", "temp_bar.png"),width = 650, height = 500)

hist(wthh$temp,
     breaks = 50,
     col = "peachpuff",
     main = "Histogram of Temperatures (Apr. - Sep.)",
     prob = TRUE) # prob = TRUE to show densities instead of frequencies
lines(density(wthh$temp, na.rm = T), lwd = 2, col = "chocolate3")
abline(v = mean(wthh$temp),col = "royalblue", lwd = 2)
abline(v = median(wthh$temp),col = "red", lwd = 2)
abline(v=quantile(wthh$temp,0.25, na.rm = T),col="darkgreen",lty=2)
abline(v=quantile(wthh$temp,0.75, na.rm = T),col="green",lty=2)

legend(x = "left", 
       bty = "n",
       c("Density plot", "Mean", "Median","Lower Quantile", "Upper Quantile"),
       col = c("chocolate3", "royalblue", "red", "darkgreen","green"),
       lwd = c(2, 2, 2))

dev.off()

png(here::here("out","wth", "rh_bars.png"),width = 650, height = 500)

hist(wthh$rh,
     col = "lightblue",
     main = "Histogram of Rel. Humidity (Apr. - Sep.)",
     prob = TRUE)
lines(density(wthh$rh, na.rm = T), lwd = 2, col = "chocolate3")
abline(v = mean(wthh$rh, na.rm = T),col = "royalblue", lwd = 2)
abline(v=quantile(wthh$rh,0.25, na.rm = T),col="darkgreen",lty=2)
abline(v=quantile(wthh$rh,0.75, na.rm = T),col="green",lty=2)
abline(v = median(wthh$rh, na.rm = T),col = "red", lwd = 2)
legend(x = "left", 
       bty = "n",
       c("Density plot", "Mean", "Median","Lower Quantile", "Upper Quantile"),
       col = c("chocolate3", "royalblue", "red", "darkgreen","green"),
       lwd = c(2, 2, 2))

dev.off()


hist(wthd$rain[wthd$rain>0], 
     breaks = 60, 
     col = "blue",
     main = "Freq. of rain > 0 mm (Apr. - Sep.)",
     xlab = "Rain (mm/hour)")


# Ridgeplots
library(ggridges)
ttl <- "Temperatures in NW Washington state per month"
wthd %>% 
  mutate(yr = year(date),
         mon = month(date)) %>% 
  dplyr::filter(yr>2000) %>%
  dplyr::filter(mon %in% c(9:12, 1:5)) %>%
  mutate(mon = factor(mon, levels = rev(c(9:12, 1:5)))) %>% 
  ggplot(  aes(x = temp, y = factor(mon, ordered = T)))+
  geom_density_ridges_gradient(
    aes(fill = ..x..), scale = 2, size = 0.2
  ) +
  scale_fill_gradientn(
    colours = c("#0D0887FF", "#CC4678FF", "#F0F921FF"),
    name = "Temp. (°C)"
  )+
  ggtitle(ttl)+
  xlab("Temperature (°C)")+
  ylab("Month")+
  theme_ridges(font_size = 13, grid = TRUE)

ggsave(here::here("out","wth", paste0(ttl, ".png")),dpi = 500)

ttl <- "Winter temperatures in NW Washington state per station"
wthd %>% 
  mutate(yr = year(date),
         mon = month(date)) %>% 
  dplyr::filter(yr>2000) %>%
  dplyr::filter(mon %in% c(11:12, 1:2)) %>%
  mutate(mon = factor(mon, levels = rev(c(9:12, 1:5)))) %>% 
  ggplot(  aes(x = temp, y = factor(stna, ordered = T)))+
  geom_density_ridges_gradient(
    aes(fill = ..x..), scale = 2, size = 0.2
  ) +
  scale_fill_gradientn(
    colours = c("#0D0887FF", "#CC4678FF", "#F0F921FF"),
    name = "Temp. (°C)"
  )+
  ggtitle(ttl)+
  xlab("Temperature (°C)")+
  ylab("Weather Station")+
  theme_ridges(font_size = 13, grid = TRUE)

ggsave(here::here("out","wth", paste0(ttl, ".png")),
       width = 8,dpi = 500)

ttl <- "Winter temperatures in NW Washington state per year"
wthd %>% 
  mutate(yr = year(date),
         mon = month(date)) %>% 
  dplyr::filter(yr>2014) %>%
  dplyr::filter(mon %in% c(11:12, 1:2)) %>%
  mutate(mon = factor(mon, levels = rev(c(9:12, 1:5)))) %>% 
  ggplot(  aes(x = temp, y = factor(yr, ordered = T)))+
  geom_density_ridges_gradient(
    aes(fill = ..x..), scale = 2, size = 0.2
  ) +
  scale_fill_gradientn(
    colours = c("#0D0887FF", "#CC4678FF", "#F0F921FF"),
    name = "Temp. (°C)"
  )+
  ggtitle(ttl)+
  xlab("Temperature (°C)")+
  ylab("Month")+
  theme_ridges(font_size = 13, grid = TRUE)

ggsave(here::here("out","wth", paste0(ttl, ".png")),dpi = 500)




# Correlation of weather
library("PerformanceAnalytics")
conflict_prefer("legend", "graphics")

load(file =  here::here("dat", "wth", "weather_day.RData"))

mtdt <- 
  read_csv( file = here("dat", "wth", "wth_loc_agwet.csv"))
 

wthd <-
  left_join(wthd, mtdt, by = "stna")

wthd$stna <-forcats::fct_reorder(wthd$stna, wthd$lat, max)
forcats::fct_rev(forcats::fct_reorder(wthd$stna, wthd$lat, max))

months_selected <- c(10:12, 1:3)
ttl <- "Correlations of all weather variables"

wth_sub <-
wthd %>% 
  mutate(yr = year(date),
         mon = month(date)) %>% 
  dplyr::filter(yr>2014) %>%
  dplyr::filter(mon %in% months_selected) %>%
  select( stna, yr, mon,date,  temp,rh,  rain, rain, solrad, wspd,soilt8in ) %>% 
  mutate(mon = factor(mon, levels = rev(months_selected))) %>%   
  pivot_longer(cols = c(temp,rh,  rain, solrad, wspd,soilt8in), 
               names_to = "vars") %>% 
  pivot_wider(names_from  = stna)





png(here::here("out","wth", "cors", paste0("corelation", ".png")),
    width = 900, height =700)

chart.Correlation(select(wth_sub, levels(wthd$stna)), histogram=FALSE, pch=19)

dev.off()

shell.exec(here::here("out","wth", paste0(ttl, ".png")))


vars <- c("temp","rh",  "rain", "solrad", "wspd","soilt8in")
for (i in vars) {
  
  png(here::here("out","wth", "cors", paste0("corelation ", i, ".png")),
      width = 900, height =700)
  wth_sub%>% 
    
    filter(vars == i) %>% 
    select(levels(wthd$stna)) %>% 
  chart.Correlation(., histogram=TRUE, pch=19)
  dev.off()
}

source("http://www.sthda.com/upload/rquery_cormat.r")
rquery.cormat(wth_sub)




