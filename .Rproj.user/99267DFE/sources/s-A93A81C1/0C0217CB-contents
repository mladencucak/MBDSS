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
(wth <- wthls[[3]])

 

# Probabilities for sporulation onset
load(here("scr/model/cu_probs.RData"))
tb
# # A tibble: 2 × 4
#     mng   means  sdev   num
#   1 high  1829.  170.     5
#   2 low   1620.  126.     7

cu <- ifelse(wth$temp > 0 & wth$temp < 7.2, 1, 0)
wth$cusum <- cumsum(cu)

wth$prob_high <- pnorm(wth$cusum,
                tb[tb$mng == "high", "means"] %>% pull,
                tb[tb$mng == "high", "sdev"] %>% pull)
wth$prob_low <- pnorm(wth$cusum,
               tb[tb$mng == "low", "means"] %>% pull, 
               tb[tb$mng == "low", "sdev"] %>% pull)

start.prob <- .01

wth$start_high <-wth[which.min(abs(wth$prob_high - start.prob)), "datetime" ]%>% pull()
wth$start_low <-wth[which.min(abs(wth$prob_low - start.prob)),"datetime" ]%>% pull()


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
temp_thresh <- 2
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
wet_dur <- ifelse(rh >= rh_thresh | rain> rain_thresh, 1,0)

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

 
 

(p1 <- 
wth %>%
  
  mutate(inf = ifelse(inf == 0 ,NA, inf)) %>%
  group_by(mng) %>% 
  mutate(inf = ifelse(datetime <spore_start ,NA, inf)) %>%  
  mutate( start_lab  = substring(as.character(as.Date(spore_start )), 6),
          # start_lab =   paste(month.abb[as.numeric(strsplit( start_lab, "-")[[1]][[1]])],
          #                     strsplit( start_lab, "-")[[1]][[2]], sep = "-")
  ) %>% 
  # mutate(col_inf = ifelse(inf < .2, "green",
  #                         ifelse(inf >=.2& inf<med_max, "orange", 
  #                                ifelse(inf >= .35, "red", "gray")))) %>% 
  filter(doy<200) %>% 
  ggplot(aes(datetime, inf))+
  geom_tile(aes(x=datetime,y=50,fill=cut(inf,3)),height=100,alpha=0.2) +
  # geom_tile(aes(x=datetime,y=50,fill=col_inf, group = 1),height=100,alpha=0.4)+
  # geom_rect(aes(xmin=spore_start, xmax=max(datetime), ymin=0, ymax=safe_max*100), 
  #           fill="#99c140") +
  # geom_rect(aes(xmin=spore_start, xmax=max(datetime), ymin=safe_max*100, ymax=med_max*100), 
  #           fill="#e7b416") +
  # geom_rect(aes(xmin=spore_start, xmax=max(datetime), ymin=med_max*100, ymax=100), 
  #           fill="#cc3232") +
  # # scale_fill_manual(values = c("#99c140",   "#e7b416", "#cc3232"))+
  
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
# load the infection model 
load( here("scr/model/inf_model.RData"))

# cumulative probability of sporulation parameters
load(here("scr/model/cu_probs.RData"))
spor.model = tb
inf.model = mod

# The model is wrapped into a function for easier apllication to lager data set
MBRisk <- function(wth, 
                   spor.model, 
                   inf.model ){
  # wth <- wthls[[1]]
  
  # Set min thresholds to enable model runs
  rh_thresh <- 90
  temp_thresh <- 2
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
  
  
  # Probabilities for sporulation onset
  # Works as a lookup table, estiamting probability of sporulation onset based on 
  
  cu <- ifelse(wth$temp > 0 & wth$temp < 7.2, 1, 0)
  wth$cusum <- cumsum(cu)
  
  wth$prob_high <- pnorm(wth$cusum,
                         spor.model[spor.model$mng == "high", "means"] %>% pull,
                         spor.model[spor.model$mng == "high", "sdev"] %>% pull)
  wth$prob_low <- pnorm(wth$cusum,
                        spor.model[spor.model$mng == "low", "means"] %>% pull, 
                        spor.model[spor.model$mng == "low", "sdev"] %>% pull)
  
  start.prob <- .01
  
  wth$start_high <-wth[which.min(abs(wth$prob_high - start.prob)), "datetime" ]%>% pull()
  wth$start_low <-wth[which.min(abs(wth$prob_low - start.prob)),"datetime" ]%>% pull()
  
  
  
  
  
  # conditions for infection
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
  wth$inf <- plogis(predict(inf.model, newdata=dff))
  
  # dff$inf <-ifelse(dff$wet_dur == 0,0, dff$inf)
  
return(wth)
  }


system.time(MBRisk(wth, tb, mod))


############################################################
# Model evaluation
############################################################

 
# Add traficlight system using the risk estimates proposed by the author of paper
# Turn the code into a function and plot all of these to get an idea if it is ok
# Calculate the average date of the initial predicted sporulation onset 
# Observed forecasted data

# Temps above 65F during infection period times 


for (i in 1:length(wthls)) {
  x <- wthls[[i]]
  # (x <- wthls[[12]])

  wthls[[i]] <-  MBRisk(x,tb, mod)

  print(paste(i,",", round(i/length(wthls),3)))
  done <- i
}

 


beepr::beep()

# check for temperatures higher than 18.33 C
(hightemps <- 
wthls %>% 
  bind_rows() %>% 
  separate(envir, into = c("stna", "season")) %>% 
  mutate(mn = month(date)) %>% 
  mutate(high_temp = ifelse(temp>18.3333, 1,0)) %>% 
  group_by(mn) %>% 
  summarise(temp = sum(high_temp)))

   
 

  # load the model 
#   load( here("scr/model/inf_model.RData"))
#   
# cl <- makeCluster(detectCores())
# clusterExport(cl, c( 
#                     "wthls",
#                     "MBRisk",
#                     "mod",
#                     "tb"
#                      ))
# 
# clusterEvalQ(cl, library("tidyverse", quietly = TRUE, verbose = FALSE))
# clusterEvalQ(cl, library("here", quietly = TRUE, verbose = FALSE))
#  
# 
# lapply(wthls, function(x){
#   MBRisk(x, mod, tb)
# })
# 
# wthlsparalel <-
#   pbapply::pblapply(wthls, function(x){
#     MBRisk(x, mod, tb)
#     }, cl = cl)
# beepr::beep()

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
    inf = sum(inf),
    start_low = unique(start_low),
    start_high = unique(start_high)
  ) %>% 
  ungroup() %>% 
  separate(envir, into = c("stna", "season"), remove = F) %>% 
  mutate(season = as.numeric(season)) %>% 
  mutate(doy = yday(date)) 
 
#Reproducble from this point 
datls <- 
dat %>%
  filter(doy <200) %>% 
  split(f= .$envir) %>% 
  lapply(., function(dff){
   df <-  dff[which.min(abs(dff$prob_low - .01)),"doy"]
   df$regime <- "low"
   df <- 
   bind_rows(df,
   data.frame(doy = dff[which.min(abs(dff$prob_high - .01)),"doy"] %>% pull(),
              regime = "high"))
   df$envir <- dff$envir[1]
   return(df)
   }) %>% bind_rows()

labs.d <- 
datls %>% 
  group_by(regime) %>% 
  summarise(`50%` = mean(doy),
            `25%` = quantile(doy, .25),
            `75%` = quantile(doy, .75),
            `5%` = quantile(doy, .05),
            `95%` = quantile(doy, .95)
  ) %>% 
  pivot_longer(cols = !regime) %>% 
  mutate(date = 
           lubridate::as_date(value, origin = "2016-01-01"),
         date.m = format(date, "%d-%m")) %>% 
  mutate( name = factor(name)) %>% 
  mutate(date.m = paste0(date.m, "(", name, ")"))
  
levels(labs.d$name) <- levels(labs.d$name)[c(2,1,3,4,5)]
 


datls %>% 
  mutate(date = 
           lubridate::as_date(doy, origin = "2016-01-01")) %>% 
  left_join(., labs.d) %>% 
  ggplot(aes(regime, date))+
  geom_boxplot(width = .5,color = "gray", notch = TRUE)+
  geom_jitter(width =.2, size = .5)+
  geom_text(aes(label = date.m, x =ifelse(regime=="high",1.5,2.5), 
                 y = date
                # color = name
                ),
            data = labs.d,
            size = 2.5,
            angle = 20
             )+
   coord_flip()+
  xlab("Management")+
  theme_bw()+
  theme(axis.title.x = element_blank())
  
ggsave(here("out/Spor_start_Dur_Intensity.png"), 
       width = 6.5, height = 3.2, dpi = 600)
shell.exec(here("out/Spor_start_Dur_Intensity.png"))


#Infection model evaluation 
datls %>% 
  group_by(regime) %>% 
  filter(doy>85) %>% head(20) %>% 
  mutate(date = 
           lubridate::as_date(doy, origin = "2016-01-01"),
         date = format(date, "%m-%d"))

 #
lsrisk <- list()
 
for (i in seq(wthls)) {
  # i = 2
  dff <- wthls[[i]]
  
  first <- which.min(abs(dff$prob_low - .01))
  duration <- 28*24
  last <- first + duration
  
  dfrisk <- 
  dff[first : last,] %>% 
    mutate(date = as_date(datetime)) %>% 
    group_by(date) %>% 
    summarise(inf = max(inf)) %>%
    mutate(risk = ifelse(inf < .2, 0,
                            ifelse(inf >=.2& inf<.35, 1, 
                                   ifelse(inf >= .35, 2, NA)))) %>% 
    group_by(risk) %>% 
    summarise(count = n())
  
  dfrisk$regime <- "low"
  
  first <- which.min(abs(dff$prob_high - .01))
  duration <- 28*24
  last <- first + duration
  
  dfriskhigh <- 
    dff[first : last,] %>% 
    mutate(date = as_date(datetime)) %>% 
    group_by(date) %>% 
    summarise(inf = max(inf)) %>%
    mutate(risk = ifelse(inf < .2, 0,
                         ifelse(inf >=.2& inf<.35, 1, 
                                ifelse(inf >= .35, 2, NA)))) %>% 
    group_by(risk) %>% 
    summarise(count = n())
  
  dfriskhigh$regime <- "high"
  lsrisk[[i]] <-  rbind(dfrisk, dfriskhigh) 
  
  }

pallet <- c("#99c140","#e7b416",  "#cc3232")

risk <- 
  lsrisk %>% 
  bind_rows() %>% 
  mutate(risk = factor(risk)) %>%
  drop_na()

labs.d <- 
  risk %>% 
  group_by(regime, risk) %>% 
  summarise(`50%` = mean(count),
            `25%` = quantile(count, .25),
            `75%` = quantile(count, .75),
            `5%` = quantile(count, .05),
            `95%` = quantile(count, .95)
  ) %>% 
  ungroup() %>% 
  pivot_longer(cols = c(`50%` ,`25%`, `75%` , `5%`, `95%`)) %>% 
  mutate(count = round(value, 1)) %>% 
  select(-value) %>% 
  mutate(date.m = paste0(count, "(", name, ")"))

# levels(labs.d$name) <- levels(labs.d$name)[c(2,1,3,4,5)]


risk %>% 
  ggplot(aes(risk, count))+
  geom_boxplot(width = .3,color = "gray", notch = TRUE)+
  geom_jitter(aes(color = risk),width =.1, size = .5)+
  # coord_flip()+
  xlab("Management")+
  scale_color_manual("Risk",
                     values = pallet, 
                     labels= c("No risk", "Medium", "High"))+
  ggrepel::geom_text_repel(aes(label = date.m, 
                x =
                  ifelse(risk ==0, 1.5, 
                         ifelse(risk == 1, 2.5,
                                ifelse(risk== 2,3.5, NA))),
                y = count,
                color = risk
  ),
  data = labs.d,
  size = 4,
  max.overlaps = 20,
  direction = "y"
  # angle = 20
  )+
  ylab("Number of days per risk category")+
  scale_x_discrete(expand = c(0,1.1))+
  theme_bw()+
  theme(axis.title.x = element_blank(),
        legend.position = "bottom",
        axis.text.x = element_blank(),
        legend.text=element_text(size=13)
        )+
  facet_wrap(~regime)

ggsave(here("out/Infection_evaluation.png"), 
       width = 10, height = 7, dpi = 600)
shell.exec(here("out/Infection_evaluation.png"))

