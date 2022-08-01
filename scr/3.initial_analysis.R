


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
    "forcats",
    "ggpubr",
    "survival"
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







#######################################################################
# Chilling units
#######################################################################
# Chill unit calculations 
load(file = here("dat/ana/chill.RData"))
#Bio data: dis_long
load(file = here::here("dat", "raw", "final", "dis_long.RData"))




# Distance between bio sites and locations  
bind_rows(ch_ls) %>% 
  separate(id, into = c("loc", "yr"), remove = FALSE) %>% 
  mutate(
    mng = ifelse(loc == "Whatcom", "low",
                 ifelse(loc == "Skagit"& yr <= 2017, "low", "high"))
  ) %>% 
  select(loc, stna, dist, yr) %>% 
  distinct() %>% 
  arrange(yr, loc)

bind_rows(ch_ls) %>% 
  separate(id, into = c("loc", "yr"), remove = FALSE) %>% 
  mutate(
    mng = ifelse(loc == "Whatcom", "low",
                 ifelse(loc == "Skagit"& yr <= 2017, "low", "high"))
  ) %>%   
select(loc, stna, dist, yr) %>% 
  distinct() %>% 
  summarise(average_distance = mean(dist),
            min_dist = min(dist),
            max_dist = max(dist))

dta <- 
bind_rows(ch_ls) %>% 
  filter(sporulated == "yes") %>%
  separate(id, into = c("loc", "yr"), remove = FALSE) %>% 
  mutate(
    mng = ifelse(loc == "Whatcom", "low",
                 ifelse(loc == "Skagit"& yr <= 2017, "low", "high"))
         )

# Mean distance between stations 
dta %>% 
  select(loc, stna, dist, yr) %>% 
  distinct() %>% 
  summarise(average_distance = mean(dist),
            min_dist = min(dist),
            max_dist = max(dist))



###########################################
# Visual exploration of cu
###########################################
 
 

# pubfig
tit <-
  "Relation of chilling units (0 - 7.2 C)  and initial sporulation date"

dta_mod <- 
dta %>% 
  group_by(mng, chill_start) %>% 
  mutate(cv =  sd(cu_) / mean(cu_) * 100,
         cv = round(cv,2),
         lab_pos = ifelse(mng == "high", max(cu_)+80, min(cu_)-80)) %>% 
  ungroup() %>% 
  group_by(chill_start) %>% 
  
  mutate(cv_overall =  sd(cu_) / mean(cu_) * 100,
         cv_overall = paste0("Overal CV =\n ",round(cv_overall,2), "(", min(cu_),", ", max(cu_),")"),
         lab_pos_overall = min(cu_)-200,
         mng = stringr::str_to_title(mng) ) %>%   
group_by(chill_start, mng) %>% 
  
  mutate(cv_group =  sd(cu_) / mean(cu_) * 100,
         cv_group = paste0(round(cv_group,2), "(", min(cu_),", ", max(cu_),")")
         )  


dta_mod %>% 
  ggplot(aes(chill_start,cu_, color = mng, group1 = chill_start, group2= mng))+
  geom_jitter(position = position_dodge(.65))+
  geom_text(aes(chill_start, lab_pos, label = cv_group, group1 = mng))+
  geom_text(aes(chill_start, lab_pos_overall, label = cv_overall), color = "black")+
  geom_boxplot( width = .4,position = position_dodge(.65), color = "darkgray")+
  geom_jitter(position = position_dodge(.65))+
  
  labs(
    x = "Start date of the accumulation period",
    y = "Chilling units (0 – 7.2 ˚C)",
    # title = tit,
    # subtitle = "Coeficient of variation between accumulations of different methods",
    color = "Management"
  ) +
  geom_vline(xintercept=seq(1.5, length(unique(dta_mod$chill_start))-0.5, 1), 
             lwd=.5, colour="gray")+
  theme_bw()+
  theme(
    legend.position = c(.9,.85),
    axis.text.x = element_text(size = 10),
    # axis text size
    axis.text.y = element_text(vjust = 0.2),
    # axis text alignment
    axis.ticks = element_line(size = 0.4),
    axis.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12)
  )

ggsave(
  filename = here::here("out", paste0(tit, ".png")),
  width = 8,
  height = 6.5,
  dpi = 300
)
shell.exec(here::here("out", paste0(tit, ".png")))

 
#######################################################
# Probability of sporulation initiation 
###########################################################

 dta_chill <-
  dta %>% 
  filter(chill_start == "11-15") %>% 
  select(cu_, mng, spor_start_prop) %>% 
  arrange(mng)

dta_chill %>% 
  ggplot(aes(x = cu_, y = spor_start_prop, col = mng)) +
  theme_bw()+geom_point()

dta_chill %>% 
  ggplot(aes(x = mng, y = cu_, col = mng)) +
  theme_bw()+geom_point()

(tb <- 
    dta_chill %>% 
    group_by(mng) %>% 
    summarise(means = mean(cu_),
              sdev = sd(cu_),
              num = n()))
# # A tibble: 2 × 4
#     mng   means  sdev   num
#   1 high  1829.  170.     5
#   2 low   1620.  126.     7


# probability of the apothecia presence is calculated for the two management levels 
# normal distribution is assumed. Such assumption and model parameters 
# can be further substantiated as the data set is growing  

save(tb, file = here("scr/model/cu_probs.RData"))


# Visualise distribution
probdf.long <- 
  probdf %>% 
  pivot_longer(cols = c(cu_high, cu_low), names_to = "mng",values_to = "calc")

p.hi <- 
dta_chill %>% 
  ggplot()+
  # geom_density(aes(cu_high,  ), data =probdf.long )+
  geom_density(aes(cu_low ), data =probdf , color = "red")+
  geom_density(aes(cu_),linetype = "dotted") +theme_bw()
p.lo <- 
dta_chill %>% 
  ggplot()+
  # geom_density(aes(cu_high,  ), data =probdf.long )+
  geom_density(aes(cu_high ), data =probdf , color = "red")+
  geom_density(aes(cu_),linetype = "dotted") +theme_bw()

ggpubr::ggarrange(plotlist = list(p.hi, p.lo)) 
 

shapiro.test(dta_chill[dta_chill$mng == "high", "cu_"] %>% pull)
shapiro.test(dta_chill[dta_chill$mng == "low", "cu_"] %>% pull)
plot(density(dta_chill[dta_chill$mng == "high", "cu_"] %>% pull))
plot(density(dta_chill[dta_chill$mng == "low", "cu_"] %>% pull))
## Plot using a qqplot
qqnorm(dta_chill[dta_chill$mng == "high", "cu_"] %>% pull)
qqline(dta_chill[dta_chill$mng == "high", "cu_"] %>% pull, col = 2)
qqnorm(dta_chill[dta_chill$mng == "high", "cu_"] %>% pull)
qqline(dta_chill[dta_chill$mng == "high", "cu_"] %>% pull, col = 2)


 
 

#######################################################
# Plot all data
#######################################################
load(file = here("dat/ana/wthd_ls.RData"))


linesize <- .35 #size of the line for variables

wthd_ls %>% bind_rows() %>% 
  mutate(rain = ifelse(rain == 0, NA, rain)) %>% 
  ggplot() +
  geom_point(aes(germ_start, 46, shape = sporulated), size = 1.5)+
  geom_point(aes(spore_start, 46),shape =2, size = 1.5)+
  geom_point(aes(spore_end, 46),shape =2, size = 1.5)+
  # geom_point(aes(date, condday), shape = 3,size = .55)+
  geom_line(
    aes(
      x = date,
      y = rh,
      colour = "Relative humidity (%)"
    ),
    size = linesize
  ) +
  geom_line(aes(
    x = date,
    y = temp,
    colour = "Temperature (˚C)",
  ),
  size = linesize) +
  geom_line(aes(
    x = date,
    y = wtemp_10 ,
    colour = "Rolling mean temp. (˚C)"
  ),
  linetype = "dashed",
  size = linesize) +
  geom_col(
    aes(date,
        rain,
        fill = "Total precipitation (mm/day)"),
    size = 1.2 ,
    inherit.aes = TRUE,
    width = 0.8
  ) +
  scale_colour_manual(
    "Daily weather:",
    values = c(
      "Relative humidity (%)" = "#0EBFE9",
      "Temperature (˚C)" = "#ED2939",
      "Rolling mean temp. (˚C)"= "darkred"
    )
  ) +
  scale_fill_manual(name= NA, values = c("Total precipitation (mm/day)" = "blue")
  ) +
  scale_y_continuous( breaks = seq(0,100,10), limits = c(-1,100),expand = c(0, 0))+
  scale_x_date(date_labels = "%b", date_breaks = "1 month",expand = c(.03, .04))+
  facet_wrap( ~ id, scales = "free", ncol = 1,strip.position="right") +
  guides(fill=guide_legend(nrow=1,byrow=TRUE,title.position = "top"),
         color=guide_legend(nrow=1,byrow=TRUE,title.position = "top"),
         linetype=guide_legend(nrow=1,byrow=TRUE,title.position = "top"))+
  theme_bw()+
  theme(
    text = element_text(size=10.8),
    # axis.text.x = element_blank(),
    # axis.ticks.x = element_blank(),
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    # axis.text.y = element_blank(),
    # axis.ticks.y = element_blank(),
    legend.position = "top",
    strip.background =element_rect(colour = "black",size=0.75), 
    strip.placement = "outside",
    legend.title = element_blank(),
    panel.grid.minor =   element_blank(),
    panel.grid.major =   element_line(colour = "lightgray",size=linesize),
    plot.margin = margin(0, 0, 0, 0, "cm"),
    panel.spacing.y=unit(0, "lines"),
    legend.box.spacing = unit(0.1, 'cm'),
    legend.margin = margin(0, 0, 0, 0, "cm"),
  )

ggsave(filename= here::here("out", "all_wth_dis.png"),
       width = 6.5, height =12, dpi = 200)

shell.exec(here("out", "all_wth_dis.png"))




####################################################
# plot months around periods
#####################################################
#Bio data: dis_long
load(file = here::here("dat", "raw", "final", "dis_long.RData"))
load(file = here("dat/ana/wthd_ls.RData"))


linesize <- .45 #size of the line for variables

top <- 30
apodot <- 3
apodot_pos <- -6
ReLab <- function(flab ,fvar ){
  paste(flab,
        substring(fvar, 6, 10))
}

# Add
lab_date_g <- 3
lab_date <- 3.7

# Extract additional data for pathogen development 
summary <- 
dis_long %>% 
  group_by(id) %>% 
  summarise(
    first_out = min(date),
    last_out = max(date)
  )
summary <- 
dis_long %>% 
  filter(stage %in% "germination") %>% 
  group_by(id) %>%
  mutate(counts_cum = cumsum(counts)) %>% 
  filter(counts_cum == max(counts_cum)) %>% 
  filter(row_number()==1) %>% 
  select(id, date) %>% 
  rename(germ_end = date) %>% 
  left_join(summary, ., by = "id")
  

(pl <- 
    wthd_ls %>% bind_rows() %>% 
    
    separate(id, c("loc", "yr"), remove = FALSE) %>% 
    mutate(doy =yday(date)) %>% 
    mutate(across(sporulation,~ ifelse(dplyr::lag(.) == 0 & . == 0 & lead(.) == 0, NA, .))) %>% 
    mutate(sporulation = ifelse(date< spore_start, NA, sporulation)) %>% 
    mutate(rain = ifelse(rain == 0, NA, rain)) %>%
    mutate(rain = ifelse(rain >top, top, rain)) %>%
    mutate(mon = month(date),
           yr = factor(year(date),
                       levels = seq(min(year(date)), max(year(date)))),
           frosthr = as.numeric(frosthr),
           frost = ifelse(frosthr >0, top-9, NA),
           spore_start = ifelse(sporulated == "yes", spore_start, NA) %>% as_date(),
           spore_end = ifelse(sporulated == "yes", spore_end, NA)%>% as_date()
    ) %>%
    left_join(., summary, by = "id") %>% 
    distinct() %>%
    dplyr::filter(mon < 5) %>%
    ggplot() +
    # Germination and sporulation areas
    geom_rect(aes(xmin = germ_start, xmax = germ_end, ymin = -Inf, ymax = Inf),
              fill=gray.colors(12)[12], alpha = 0.9)+
    geom_rect(aes(xmin = spore_start, xmax = spore_end, ymin = -Inf, ymax = Inf),
              fill=gray.colors(12)[10], alpha = 0.9)+
    geom_vline(aes(xintercept = germ_start), color = "black", size  = .2)+
    geom_vline(aes(xintercept = spore_start), color = "black", size  = .2)+
    geom_vline(aes(xintercept = spore_end), color = "black", size  = .2)+
    geom_vline(aes(xintercept = last_out), color = "black", size  = .4, linetype = "dotted")+
    
    # Boost the zero line
    geom_hline(aes(yintercept = 0), color =gray.colors(12)[4] )+
    
    # Environmental variables
    geom_line(aes( x = date, y = temp,colour = "Mean daily temp.(˚C)"),
              size = linesize) +
    geom_line(aes(x = date, y = wmintemp_10,colour = "10-day rolling mean daily min. temp.(˚C)"),
              size = linesize, linetype = "twodash") +
    geom_col(aes(date,rain,fill = "Total daily precipitation (mm/day)"), 
             size = .9 ,inherit.aes = TRUE, width = 0.4) +
    # Negative temperatures
    geom_point(aes(date, frost, group = frosthr, size = frosthr), shape = 20)+
    
    # Labels for germination and sporulation dates    
    geom_text(aes(germ_start - lab_date_g, apodot_pos, label = ReLab( "G",germ_start)), size = apodot,show.legend = FALSE)+
    geom_text(aes(spore_start - lab_date, apodot_pos, 
                  label = ReLab(  ifelse(spore_start == spore_end, "AS-E", "AS"),spore_start)),size = apodot,show.legend = FALSE)+
    geom_text(aes(spore_end + lab_date, apodot_pos, 
                  label = ifelse(spore_start == spore_end, NA, ReLab( "AE", spore_end))),size = apodot,show.legend = FALSE)+
    geom_text(aes(germ_start+ 30, apodot_pos, 
                  label = ifelse(is.na(spore_start),"No apothecia presence", NA)),size = apodot,show.legend = FALSE)+
    # Sporulation percentages
    geom_point(aes(date, ifelse(sporulation > top, top-6, sporulation) , 
                   group = sporulation, colour = "Sporulating apothecia (%)"), shape =17, size = 2, alpha = .8) +
    geom_text(aes(date-1, top-3, 
                  label =ifelse(!is.na(sporulation), paste0(sporulation, "%"), NA),
                  colour = "Apothecia percentage"
    ),size = apodot,show.legend = FALSE)+
    #Theme and colors 
    scale_colour_manual(
      "Daily weather:",
      values = c(
        "Mean daily temp.(˚C)" = "#ED2939",
        "10-day rolling mean daily min. temp.(˚C)" = "darkred",
        "Sporulating apothecia (%)" = "darkgreen"
      )
    ) +
    scale_fill_manual(name= NA, values = c("Total daily precipitation (mm/day)" = "blue")
    ) +
    scale_y_continuous( breaks = seq(-10,top,10), limits = c(-10,top),expand = c(0, 0))+
    scale_x_date(date_labels = "%b", date_breaks = "1 month",expand = c(.03, .04))+
    
    guides(fill=guide_legend(nrow=1,byrow=F,title.position = "top",override.aes = list(linetype = 0)),
           color=guide_legend(nrow=2,byrow=F,title.position = "top",
                              override.aes = list(shape = NA,linetype ="solid" )),
           linetype=guide_legend(nrow=1,byrow=T,title.position = "top",override.aes = list(linetype ="twodash"))
           )+
    theme_bw()+
    theme(
      text = element_text(size=10.8),
      # axis.text.x = element_blank(),
      # axis.ticks.x = element_blank(),
      axis.title.y = element_blank(),
      axis.title.x = element_blank(),
      # axis.text.y = element_blank(),
      # axis.ticks.y = element_blank(),
      legend.position = "top",
      strip.background =element_rect(colour = "black",size=0.75), 
      strip.placement = "outside",
      legend.title = element_blank(),
      panel.grid.minor =   element_blank(),
      panel.grid.major =   element_line(colour = "lightgray",size=linesize),
      plot.margin = margin(0, 0, 0, 0, "cm"),
      panel.spacing.y=unit(0, "lines"),
      legend.box.spacing = unit(0.1, 'cm'),
      legend.margin = margin(0, 0, 0, 0, "cm"),
    )+
    facet_wrap(yr ~ loc, scales = "free", ncol = 1,strip.position="right") 
  )


pl+  
  facet_wrap(yr ~ loc, scales = "free", ncol = 1,strip.position="right") 
ggsave(filename= here::here("out", "all_wth_dis_zoom_stations.png"),
       width = 10, height =15, dpi = 450)
shell.exec(here("out", "all_wth_dis_zoom_stations.png"))


pl+  
  facet_wrap(loc ~ yr, scales = "free", ncol = 1,strip.position="right") 
ggsave(filename= here::here("out", "all_wth_dis_zoom_stations_years.png"),
       width = 10, height =15, dpi = 350)


shell.exec(here("out", "all_wth_dis_zoom_stations_years.png"))
shell.exec(here("out"))

  