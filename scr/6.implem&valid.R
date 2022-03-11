#######################################################################
# Model implementation
#######################################################################

# Inputs: 
# Variables used temp(C), rh(%), rain(mm/l)
# Temporal resolution: hourly 

# Tow models 
# 
# The accumulation (cumulative sum) of chill hours begins on 15 of November previous year

# Chill hours are are thermal units: 0-7.2 C 
# The minim number of chill hours will indicate the beginning of sporulation 
# period in two distinct management levels
# Low management = typically starts earlier 
# High management - starts later 

# Probabilities of sporulation occurrence are then calculated based on 
# the distribution of historically observed sums of chill hour units 


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
    "tidyr",
    "pbapply",
    "parallel",
    "zoo"
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
# Weather data 
#####################################################


load(here("dat", "wth", "weather_hourly.RData"))
# wthh$yr <- year(wthh$datetime)
# wthh$mon <- lubridate::month(wthh$datetime)
wthh$date <- as.Date(wthh$datetime)

wthh %>% 
  group_by(stna) %>% 
  summarise(start = min(date))

head(wthh)

wthh %>% 
  mutate(
    yr = year(date)
  )



leap_start <- 320
non_leap_start <- 319
biofix <- "11-15"
termination_date <- "05-01"

# Create environments for each season from 15th of November 
wthh <-
  wthh %>% 
  select(stna, datetime, temp, rh, rain, date) %>% 
  mutate(
    yr = year(datetime),
    #Find 15th November for each year
    start_accumulation_day = ifelse(leap_year(year(date)),
                                    leap_start, 
                                    non_leap_start),
    season = ifelse(yday(datetime) >= start_accumulation_day, yr+1, yr )
  ) %>% 
  select(-start_accumulation_day) %>% 
  unite("envir", c("stna", "season"), remove = FALSE)  


# Make a list of environments  
wthls <- split(wthh, wthh$envir)


#Clean up environments 
# Remove each environment which does not have the satisfactory data set
length(wthls) #check at the beginning

names(wthls)

for (i in names(wthls) ) {
  # i <- names(wthls)[1]
  x <-
    wthls[[i]]  
  
  
  
  x$doy <-lubridate::yday(x$datetime)
  x
  
  
  if(
    #Remove each envir that does not contain biofix or end date 
    !leap_start %in% x$doy | !non_leap_start %in% x$doy |
    x$doy[1] > leap_start | x$doy[1] > non_leap_start |
    #Remove each envir that does not have two different seasons
    length(unique(x$yr))!=2 | 
    # If the season does not contain envir_duration number of days
    length(unique(x$date)) < 140
  ){
    wthls[[i]] <- NULL
  } else {
    start_date <- as.Date(paste(unique(x$yr)[1], biofix, sep = "-"))
    end_date <-
      as.Date(paste(unique(x$yr)[2], termination_date, sep = "-"))
    x <- 
      filter(x, date >= start_date & date <= end_date) %>% 
      select(-stna, -doy, -season, -yr)
    wthls[[i]] <- x
  }
}

# Remove environments with missing data
na_prop <- function(vec){
  sum(is.na(vec))/length(vec)
}
for (i in names(wthls) ) {
  x <-wthls[[i]]  
  if(na_prop(x$temp)> 0 | na_prop(x$rh)> 0)wthls[[i]] <- NULL
}

length(wthls)



##################################################
# Model implementation
# 1. Sporulation conditions initiation
################################################
# Sample weather data
wth <- wthls[[3]]

 
# Probabilities for sporulation onset
# Works as a lookup table, estiamting probability of sporulation onset based on 
# nearest chilling unit vallue (cu_)
(load(here("scr/model/cu_probs.RData")))

dtlow = data.table(probdf, val = probdf$cu_low ) 
setattr(dtlow, "sorted", "cu_low")  # let data.table know that variable is sorted
setkey(dtlow, val) # sorts the data
dthigh =  data.table(probdf, val = probdf$cu_high  ) 
setattr(dthigh, "sorted", "cu_high")   
setkey(dthigh, val)  


wth$cu <-  ifelse(wth$temp > 0 & wth$temp < 7.2, 1, 0)
wth$cusum <- cumsum(wth$cu)
wth$cu <-NULL

# if(max(x$cusum)<min(probdf$cu_low) )
wth$prob_low <- 0
wth[wth$cusum < min(probdf$cu_low) ,"prob_low" ] <- 0
wth[ wth$cusum > max(probdf$cu_low) ,"prob_low" ] <- 1

wthx <- wth[wth$cusum >= min(probdf$cu_low) & wth$cusum <= max(probdf$cu_low) , ]

wthx$prob_low  <- 
  sapply(wthx$cusum, function(y){dtlow[J(y), roll = "nearest"][ 1,1] %>% pull()})

wth[wth$cusum >= min(probdf$cu_low) & wth$cusum <= max(probdf$cu_low), "prob_low"] <- 
  wthx$prob_low



wth$prob_high <- 0
wth[wth$cusum < min(probdf$cu_high) ,"prob_high" ] <- 0
wth[ wth$cusum > max(probdf$cu_high) ,"prob_high" ] <- 1

wthx <- wth[wth$cusum >= min(probdf$cu_high)& wth$cusum <= max(probdf$cu_high) , ]

wthx$prob_high  <- 
  sapply(wthx$cusum, function(y){dthigh[J(y), roll = "nearest"][ 1,1] %>% pull()})

wth[wth$cusum >= min(probdf$cu_high)& wth$cusum <= max(probdf$cu_high) ,"prob_high" ] <- 
  wthx$prob_high

wth$start_low <- wth[wth$cusum >= min(probdf$cu_low),"datetime"][1,] %>% pull()
wth$start_high <- wth[wth$cusum >= min(probdf$cu_high),"datetime"][1,] %>% pull()
wth[3000:3200, ] %>% 
  mutate(dte = as.Date(datetime)) %>% 
  select(dte, cusum, prob_high, prob_low) %>% 
  group_by(dte) %>% 
  summarise_all(mean)

##################################################
# Model implementation
# 1. Sporulation conditions initiation
# 2. Infection Model
##########################################################
# Starts when minimum conditions for the minimum conditions for 
# the sporulation have been met, set to  - prob >.01  
  

# load the model 
load( here("scr/model/inf_model.RData"))

# Set min thresholds to enable model runs
rh_thresh <- 90
temp_thresh <- 6
rain_thresh <- .2

# Extract variables as vectors for speed
wth[["rain"]] -> rain
if ("rhum" %in% names(wth))  wth[["rhum"]] -> rh
if ("rh" %in% names(wth))  wth[["rh"]] -> rh
wth[["temp"]] -> temp


# This function to infill missing values to let the model run
infill_gap <- 12

if (sum(is.na(with(wth, rain, temp, rhum))) > 0) {
  temp <-
    round(zoo::na.spline(temp, na.rm = FALSE, maxgap = infill_gap), 1)
  rh <-
    round(zoo::na.spline(rh, na.rm = FALSE, maxgap = infill_gap), 0)
  rh  <- sapply(rh, function(x) ifelse(x > 100, x <- 100, x))
  # Rain is infilled based on min rh wethens threshold to let the model run
  rain <- ifelse(rh>=rh_thresh, .2 )
}

if (sum(is.na(with(wth, rain, temp, rhum))) > 0) {
  stop(print("The sum of NAs is more than 7! Check your weather data."))
}



# conditions for sporulation
wet_dur <- ifelse(rh >= 90 | rain> rain_thresh, 1,0)

criteria <- as.numeric(temp >= temp_thresh & wet_dur == 1)


# criteria  <- c(0,1,0,1,0,0,1,1,0,0,0,1,1,0,0,0)
#The accumulation breaks if the conditions aren't met for more than infstop hours
infstop <- 2 + 1
criteria <- c(criteria, rep(0, infstop))

for (k in c(1:c(length(criteria)-infstop))){
  # k = 2
  if(criteria[k] == 1& criteria[k + infstop] == 1 ) criteria[k : c(k + infstop)] <- 1
}
criteria <-  criteria[1:c(length(criteria)-infstop)]
#

# cumulative sum of hours that meet the criteria for sporulation with restart at zero
(criteria_sum <-
    stats::ave(criteria, cumsum(criteria == 0), FUN = cumsum)
)
dff <- data.frame(temp = temp, wet_dur=criteria_sum)
wth$inf <- plogis(predict(mod, newdata=dff))

# dff$inf <-ifelse(dff$wet_dur == 0,0, dff$inf)





wth <-
  mutate(wth,doy = yday(datetime)) %>% 
  filter(doy <200) %>% 
  pivot_longer(cols = c(prob_low, prob_high),
               names_to = "mng" ) %>% 
  mutate(mng=ifelse(mng == "prob_low", "Low", "High"),
         inf = inf * 100) %>% 
  pivot_longer(cols = c(start_low, start_high),
               values_to = "spore_start",
               names_to = "spore_start_mng" )  %>% 
  mutate(spore_start_mng=ifelse(spore_start_mng == "start_low", "Low", "High")) %>% 
         
  mutate(value = 100 *value)


#Parameteres
safe_max <- .2
med_max <- .35
text_size <- 14
 
 
format.mmdd <- function(x, format = "%b-%d", ...) format(as.Date(x), format = format, ...)


wth %>% 
  mutate(inf = ifelse(inf == 0 ,NA, inf)) %>%
  group_by(mng) %>% 
  mutate(inf = ifelse(datetime <spore_start ,NA, inf)) %>%  
  mutate( start_lab  = substring(as.character(as.Date(spore_start )), 6),
          # start_lab =   paste(month.abb[as.numeric(strsplit( start_lab, "-")[[1]][[1]])],
          #                     strsplit( start_lab, "-")[[1]][[2]], sep = "-")
  ) %>% 
  mutate(col_inf = ifelse(inf < 2, "green",
                          ifelse(inf >=2& inf<.4, "orange", 
                                 ifelse(inf >= .4, "red", "gray")))) %>% 
  filter(doy<200) 
 
 

(p1 <- 
wth %>%
  
  mutate(inf = ifelse(inf == 0 ,NA, inf)) %>%
  group_by(mng) %>% 
  mutate(inf = ifelse(datetime <spore_start ,NA, inf)) %>%  
  mutate( start_lab  = substring(as.character(as.Date(spore_start )), 6),
          # start_lab =   paste(month.abb[as.numeric(strsplit( start_lab, "-")[[1]][[1]])],
          #                     strsplit( start_lab, "-")[[1]][[2]], sep = "-")
  ) %>% 
  mutate(col_inf = ifelse(inf < .2, "green",
                          ifelse(inf >=.2& inf<.4, "orange", 
                                 ifelse(inf >= 40, "red", "gray")))) %>% 
  filter(doy<200) %>% 
  ggplot(aes(datetime, inf))+
  # geom_tile(aes(x=datetime,y=50,fill=cut(inf,3)),height=100,alpha=0.2) +
  # geom_tile(aes(x=datetime,y=50,fill=col_inf, group = 1),height=100,alpha=0.4)+
  geom_rect(aes(xmin=spore_start, xmax=max(datetime), ymin=0, ymax=safe_max*100), 
            fill="#99c140") +
  geom_rect(aes(xmin=spore_start, xmax=max(datetime), ymin=safe_max*100, ymax=med_max*100), 
            fill="#e7b416") +
  geom_rect(aes(xmin=spore_start, xmax=max(datetime), ymin=med_max*100, ymax=100), 
            fill="#cc3232") +
  # scale_fill_manual(values = c("#99c140",   "#e7b416", "#cc3232"))+
  
  geom_line(aes(datetime, inf))+
  scale_y_continuous(limits = c(0,100))+
  geom_point(aes(spore_start, .04, color = spore_start_mng, group  =2), shape = 25, size  = 2)+
  geom_text(aes(label = start_lab , x = spore_start, 
                y = ifelse(spore_start_mng == "Low", 10,30) ),
            check_overlap = TRUE,angle = 90, size = 4.5)+
  ylab("Infection Risk(%)")+
  
  theme_bw()+
  theme(axis.title.x = element_blank(),
        legend.position = "none",
        text = element_text(size = text_size))
  )


(p2 <-
    wth %>% 
    filter(doy<200) %>% 
ggplot() +
  geom_line(aes(datetime , value, color = mng) )+
  scale_y_continuous(limits = c(0, 100))+
  geom_point(aes(  spore_start, .04, fill = spore_start_mng), 
             shape = 25, size  = 2, guide_legend = FALSE)+
    ylab("Sporulation initiated(%)")+
    scale_color_discrete(name ="Management level:" )+
    scale_fill_discrete(name ="Management level:" )+
    theme_bw()+
    theme(axis.title.x = element_blank(),
          legend.position = "bottom",
          text = element_text(size = text_size) 
          )+
    theme(panel.grid.major.y = element_line(color = "black",
                                            size = 0.5,
                                            linetype = 2))
  )

 
plotf <- egg::ggarrange(plots = list(p1, p2))

ggsave(plot = plotf,
       filename = here("scr/model/visualisation.png"),
       width = 11,
       height = 7,
       dpi = 600)
shell.exec(here("scr/model/visualisation.png"))

 
    wth %>%
      filter(doy %in% c(40:50)) %>% 

        mutate(inf = ifelse(inf == 0 ,NA, inf)) %>%
        group_by(mng) %>% 
        mutate(inf = ifelse(datetime <spore_start ,NA, inf)) %>%  
        mutate( start_lab  = substring(as.character(as.Date(spore_start )), 6),
                # start_lab =   paste(month.abb[as.numeric(strsplit( start_lab, "-")[[1]][[1]])],
                #                     strsplit( start_lab, "-")[[1]][[2]], sep = "-")
        ) %>% 
        mutate(col_inf = ifelse(inf < .2, "green",
                                ifelse(inf >=.2& inf<.4, "orange", 
                                       ifelse(inf >= 40, "red", "gray")))) %>% 
        filter(doy<200) %>% 
        ggplot(aes(datetime, inf))+
        # geom_tile(aes(x=datetime,y=50,fill=cut(inf,3)),height=100,alpha=0.2) +
      geom_tile(aes(x=datetime,y=50,fill=col_inf, group = 1),height=100,alpha=0.4)+
      # geom_rect(aes(xmin=spore_start, xmax=max(datetime), ymin=0, ymax=safe_max*100), 
      #           fill="#99c140") +
      # geom_rect(aes(xmin=spore_start, xmax=max(datetime), ymin=safe_max*100, ymax=med_max*100), 
      #           fill="#e7b416") +
      # geom_rect(aes(xmin=spore_start, xmax=max(datetime), ymin=med_max*100, ymax=100), 
      #           fill="#cc3232") +
      scale_fill_manual(values = c("#99c140",   "#e7b416", "#cc3232"))+
        
        geom_line(aes(datetime, inf))+
        scale_y_continuous(limits = c(0,100))+
        geom_point(aes(spore_start, .04, color = spore_start_mng, group  =2), shape = 25, size  = 2)+
        geom_text(aes(label = start_lab , x = spore_start, 
                      y = ifelse(spore_start_mng == "Low", 10,30),
                      color = spore_start_mng),
                  check_overlap = TRUE,angle = 90, size = 4.5)+
        ylab("Infection Risk(%)")+
        
        theme_bw()+
        theme(axis.title.x = element_blank(),
              legend.position = "none")
 
    
    wth %>% 
      mutate(inf = ifelse(inf == 0 ,NA, inf)) %>%
      group_by(mng) %>% 
      mutate(inf = ifelse(datetime <spore_start ,NA, inf)) %>%  
      mutate( start_lab  = substring(as.character(as.Date(spore_start )), 6),
              # start_lab =   paste(month.abb[as.numeric(strsplit( start_lab, "-")[[1]][[1]])],
              #                     strsplit( start_lab, "-")[[1]][[2]], sep = "-")
      ) %>% 
      mutate(col_inf = ifelse(inf < 20, "green",
                              ifelse(inf >=20& inf<40, "orange", 
                                     ifelse(inf >= 40, "red", "gray")))) %>% 
      filter(doy<200) %>% 
      ggplot(aes(datetime, inf))+
      # geom_tile(aes(x=datetime,y=50,fill=cut(inf,3)),height=100,alpha=0.2) +
      geom_tile(aes(x=datetime,y=50,fill=col_inf),height=100,alpha=0.4)+
      scale_fill_manual(values = c("green", "orange", "red", "grey", "black"))+
      
      geom_line(aes(datetime, inf))+
      scale_y_continuous(limits = c(0,100))+
      geom_point(aes(  spore_start, .04, fill = spore_start_mng), shape = 25, size  = 2)+
      geom_text(aes(label = start_lab , x = spore_start, 
                    y = ifelse(spore_start_mng == "Low", 10,30),
                    color = spore_start_mng),
                check_overlap = TRUE,angle = 90, size = 4.5)+
      ylab("Infection Risk(%)")+
      
      theme_bw()+
      theme(axis.title.x = element_blank(),
            legend.position = "none")
    
    
    
ggplot(wth, aes(datetime , cusum))+
  geom_line()
# ggplot(wth, aes(datetime , prob_low))+
#   geom_line()
# 

 




##########################################################
# The model function
##########################################################

# The model is wrapped into a function for easier apllication to lager data set
MBRisk <- function(wth){
  
  
  
  dtlow = data.table(probdf, val = probdf$cu_low ) 
  setattr(dtlow, "sorted", "cu_low")  # let data.table know that variable is sorted
  setkey(dtlow, val) # sorts the data
  dthigh =  data.table(probdf, val = probdf$cu_high  ) 
  setattr(dthigh, "sorted", "cu_high")   
  setkey(dthigh, val)  
  
  
  wth$cu <-  ifelse(wth$temp > 0 & wth$temp < 7.2, 1, 0)
  wth$cusum <- cumsum(wth$cu)
  wth$cu <-NULL
  
  # if(max(x$cusum)<min(probdf$cu_low) )
  wth$prob_low <- 0
  wth[wth$cusum < min(probdf$cu_low) ,"prob_low" ] <- 0
  wth[ wth$cusum > max(probdf$cu_low) ,"prob_low" ] <- 1
  
  wthx <- wth[wth$cusum >= min(probdf$cu_low) & wth$cusum <= max(probdf$cu_low) , ]
  
  wthx$prob_low  <- 
    sapply(wthx$cusum, function(y){dtlow[J(y), roll = "nearest"][ 1,1] %>% pull()})
  
  wth[wth$cusum >= min(probdf$cu_low) & wth$cusum <= max(probdf$cu_low), "prob_low"] <- 
    wthx$prob_low
  
  
  
  wth$prob_high <- 0
  wth[wth$cusum < min(probdf$cu_high) ,"prob_high" ] <- 0
  wth[ wth$cusum > max(probdf$cu_high) ,"prob_high" ] <- 1
  
  wthx <- wth[wth$cusum >= min(probdf$cu_high)& wth$cusum <= max(probdf$cu_high) , ]
  
  wthx$prob_high  <- 
    sapply(wthx$cusum, function(y){dthigh[J(y), roll = "nearest"][ 1,1] %>% pull()})
  
  wth[wth$cusum >= min(probdf$cu_high)& wth$cusum <= max(probdf$cu_high) ,"prob_high" ] <- 
    wthx$prob_high
  
  wth$start_low <- wth[wth$cusum >= min(probdf$cu_low),"datetime"][1,] %>% pull()
  wth$start_high <- wth[wth$cusum >= min(probdf$cu_high),"datetime"][1,] %>% pull()
  
  
  
  # Set min thresholds to enable model runs
  rh_thresh <- 90
  temp_thresh <- 6
  rain_thresh <- .2
  
  # Extract variables as vectors for speed
  wth[["rain"]] -> rain
  if ("rhum" %in% names(wth))  wth[["rhum"]] -> rh
  if ("rh" %in% names(wth))  wth[["rh"]] -> rh
  wth[["temp"]] -> temp
  
  
  # This function to infill missing values to let the model run
  infill_gap <- 12
  
  if (sum(is.na(with(wth, rain, temp, rhum))) > 0) {
    temp <-
      round(zoo::na.spline(temp, na.rm = FALSE, maxgap = infill_gap), 1)
    rh <-
      round(zoo::na.spline(rh, na.rm = FALSE, maxgap = infill_gap), 0)
    rh  <- sapply(rh, function(x) ifelse(x > 100, x <- 100, x))
    # Rain is infilled based on min rh wethens threshold to let the model run
    rain <- ifelse(rh>=rh_thresh, .2 )
  }
  
  if (sum(is.na(with(wth, rain, temp, rhum))) > 0) {
    stop(print("The sum of NAs is more than 7! Check your weather data."))
  }
  
  
  
  # conditions for sporulation
  wet_dur <- ifelse(rh >= 90 | rain> rain_thresh, 1,0)
  
  criteria <- as.numeric(temp >= temp_thresh & wet_dur == 1)
  
  
  # criteria  <- c(0,1,0,1,0,0,1,1,0,0,0,1,1,0,0,0)
  #The accumulation breaks if the conditions aren't met for more than infstop hours
  infstop <- 2 + 1
  criteria <- c(criteria, rep(0, infstop))
  
  for (k in c(1:c(length(criteria)-infstop))){
    # k = 2
    if(criteria[k] == 1& criteria[k + infstop] == 1 ) criteria[k : c(k + infstop)] <- 1
  }
  criteria <-  criteria[1:c(length(criteria)-infstop)]
  #
  
  # cumulative sum of hours that meet the criteria for sporulation with restart at zero
  (criteria_sum <-
      stats::ave(criteria, cumsum(criteria == 0), FUN = cumsum)
  )
  dff <- data.frame(temp = temp, wet_dur=criteria_sum)
  wth$inf <- plogis(predict(mod, newdata=dff))
  
  # dff$inf <-ifelse(dff$wet_dur == 0,0, dff$inf)
  
  return(dff)
}


MBRisk(wth)


############################################################
# Model validation
############################################################

 
# Add traficlight system using the risk estimates proposed by the author of paper
# Turn the code into a function and plot all of these to get an idea if it is ok
# Calculate the average date of the initial predicted sporulation onset 
# Plot this on the map? 
# Observed forecasted data

# Temps above 65F during infection period times 


# for (i in 1:length(wthls)) {
#   x <- wthls[[i]]
#   # (x <- wthls[[12]])
# 
#   wthls[[i]] <-  MBRisk(wth)
# 
#   print(paste(i,",", round(i/length(wthls),3)))
#   done <- i
# }

beepr::beep()

hightemps <- 
wthls %>% 
  bind_rows() %>% 
  separate(envir, into = c("stna", "season")) %>% 
  mutate(mn = month(date)) %>% 
  mutate(high_temp = ifelse(temp>18.3333, 1,0)) %>% 
  group_by(mn) %>% 
  summarise(temp = sum(high_temp))

  mutate(hightemps, mn = factor(mn, levels =c(11, 12, 1:5)))
  
 

  # load the model 
  load( here("scr/model/inf_model.RData"))
  
cl <- makeCluster(detectCores())
clusterExport(cl, c("probdf",
                    "wthls",
                    "MBRisk",
                    "mod"
                     ))

clusterEvalQ(cl, library("tidyverse", quietly = TRUE, verbose = FALSE))
clusterEvalQ(cl, library("here", quietly = TRUE, verbose = FALSE))
clusterEvalQ(cl, library("data.table", quietly = TRUE, verbose = FALSE))


 

wthls <-
  pbapply::pblapply(wthls, MBRisk, cl = cl)
beepr::beep()

# save(wthls, file = here("out/val/wthls.RData"))
load(file = here("out/val/wthls.RData"))


dat <- 
wthls %>% 
  bind_rows() %>% 
  mutate(date = as.Date(datetime)) %>% 
  group_by(envir, date) %>% 
  summarise( 
    prob_low = mean(prob_low),
    prob_high = mean(prob_high),
    inf = sum(inf)
  ) %>% 
  ungroup() %>% 
  separate(envir, into = c("stna", "season")) %>% 
  
  mutate(doy = yday(date)) 

 dat %>% 
   pivot_longer(cols = c(prob_low, prob_high),
                names_to = "mng" ) %>% 
   mutate(mng=ifelse(mng == "prob_low", "Low Mng", "High Mng."),
          inf = inf * 100) %>% 
   
   mutate(value = 100 *value) %>% 
  ggplot()+
  geom_line(aes(doy, value), size = .01)+
  facet_wrap(~stna)


 dat %>%
   filter(doy <200) %>% 
   group_by(stna, doy) %>% 
   summarise(inf = mean(inf)) %>% 
   group_by(stna) %>% 
   mutate(cumulative_infection_risk = cumsum(inf)) %>% 
   ggplot(aes(doy, cumulative_infection_risk, color = stna)) +
   geom_line( size = .01)
   facet_wrap( ~ stna)
 
 dat %>%
   filter(doy <200) %>% 
   group_by(stna, doy) %>% 
   summarise(inf = sum(inf))
 
 
 
 