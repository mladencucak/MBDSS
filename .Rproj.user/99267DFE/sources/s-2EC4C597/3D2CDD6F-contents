###################################################
#Libraries
#####################################################

list.of.packages <-
  c(
    "dplyr",
    "ggplot2",
    "readr",
    "forcats",
    "tidyr",
    "lubridate",
    "here",
    "stringr",
    "conflicted",
    "lubridate",
    "ggpubr",
    "scales",
    "survival",
    "hnp"
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
# Statistical analysis of ascospore presence patterns 
#######################################################################

# Chill unit calculations 
load(file = here("dat/ana/chill.RData"))

dta <- 
  bind_rows(ch_ls) %>% 
  separate(id, into = c("loc", "yr"), remove = FALSE) %>% 
  mutate(
    mng = ifelse(loc == "Whatcom", "low",
                 ifelse(loc == "Skagit"& yr <= 2017, "low", "high"))
  )
spor_dur <- 
  dta %>%     
  group_by(id) %>%
  filter(chill_start == "11-15") %>%
  select(id, cu_, mng,spor_dur) %>% distinct()


#Bio data: dis_long
load(file = here::here("dat", "raw", "final", "dis_long.RData"))


dtasum <- 
  dis_long %>%
  group_by(id) %>%
  filter(stage == "sporulation") %>% 
  mutate(times = yday(date) - min(yday(date))) %>% 
  mutate(spor = sum(diff(times)*zoo::rollmean(prop,2))) %>% 
  select(id,  spor) %>%
  distinct() %>% 
  bind_rows(., tibble(id = "Skagit_2019", spor = 0))# There were no readings for sporulation for skagit  for that year


# Find the date of the first observed fully developed apothecia 
dtasum <- 
  dis_long %>%
  group_by(id) %>%
  filter(counts>0) %>% 
  filter(stage == "sporulation") %>% 
  mutate(spor_start = min(yday(date))) %>% 
  select(id, spor_start) %>%
  distinct() %>% 
  left_join(dtasum,.,by = "id") %>% 
  replace_na(list(spor_start =0))# years without sporulation will be returned as NA after joining 



dtasum <-
  left_join(spor_dur,dtasum,by = "id" ) %>% 
  separate(id, into = c( "Location","Year")) %>% 
  rename(Management = mng)

dtasum <-
  dtasum %>% 
  mutate_if(is.character, factor) %>% 
  mutate(Management = fct_rev(Management))

dtasum %>%
  filter(spor_start > 0) %>%
  select(spor_start, cu_, Management) %>% 
  distinct()
 
################################################
# Analysis 
################################################

## Start of sporulation
## Gamma GLM
l.start <-
  dtasum %>%
  filter(spor_start > 0) %>% 
  lm(spor_start ~ Management + Year + Location, data = .)
anova(l.start)
m.start <- 
   dtasum %>%
  filter(spor_start > 0) %>% 
  coxph(Surv(spor_start, event = rep(1, nrow(.) )) ~
                   Management + Year + Location,
        data= .
                 )
anova(m.start)



## Duration 
m.dur <- 
  dtasum %>%
  filter(spor_start > 0) %>%
  
lm(spor_dur  ~ Management + Year + Location, data =  .)
anova(m.dur) 
summary(m.dur)

## Cox proportional hazards model for time until event
m.dur <- 
  dtasum %>%
  filter(spor_start > 0) %>%
  coxph(Surv(spor_dur, event = rep(1, nrow(.))) ~
                 Management + Year + Location, data = .)
anova(m.dur)



# Sporulation (area under the curve)
m.spor <- lm(spor ~ Management + Year + Location, data =  dtasum)
anova(m.spor)
summary(m.spor)

m.spor <- glm(spor + .01 ~ Year + Management + Location,
              family = Gamma(link = log),
              data = dtasum)
drop1(m.spor, test = "Chisq")
scales::pvalue(as.numeric(drop1(m.spor, test = "Chisq")[[4]][2]))



###############################################
# Visualize data
###############################################

label_size = 5
year_x_lab = 6
loc_x_lab = 5



(p1<-
    dtasum %>%
    filter(spor_start > 0) %>%
    ggplot(aes(x = Management, y = spor_start)) +
    theme_bw() +
    geom_boxplot(width= .2,alpha = .3, color = "gray") +
    geom_jitter(height = 0, width = .1)+
    geom_text(aes(2.5, 65),
              label=paste0("italic(p)",scales::pvalue(as.numeric(anova(m.start)[[4]][2]))),
              parse=TRUE,
              hjust=1, size=label_size)+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())+
    ylab("Start of Sporulation(DOY)"))


(p2<-
    dtasum %>%
    filter(spor_start > 0) %>%
    ggplot(aes(x = Year, y = spor_start)) +
    theme_bw() +
    geom_boxplot(width= .2,alpha = .3, color = "gray") +
    geom_text(aes(year_x_lab, 65),
              label=paste0("italic(p)==",scales::pvalue(as.numeric(anova(m.start)[[4]][3]))),
              parse=TRUE,
              hjust=1, size=label_size)+
    geom_jitter(height = 0, width = .1)+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank()))


(p3<-
    dtasum %>%
    filter(spor_start > 0) %>%
    ggplot(aes(x = Location, y = spor_start)) +
    theme_bw() +
    geom_boxplot(width= .2,alpha = .3, color = "gray") +
    geom_text(aes(loc_x_lab, 60),
              label=paste0("italic(p)==",scales::pvalue(as.numeric(anova(m.start)[[4]][4]))),
              parse=TRUE,
              hjust=1, size=label_size)+
    geom_jitter(height = 0, width = .1)+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank()))

(p21<-
    dtasum %>%
    filter(spor_start > 0) %>%
    ggplot(aes(x = Management, y = spor_dur)) +
    theme_bw() +
    geom_boxplot(width= .3,alpha = .3, color = "gray") +
    geom_text(aes(2.5, -2),
              label=paste0("italic(p)==",scales::pvalue(as.numeric(anova(m.dur)[[4]][2]))),
              parse=TRUE,
              hjust=1, size=label_size)+
    geom_jitter(height = 0, width = .1)+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())+
    ylab("Spor. duration(days)"))

(p22<-
    dtasum %>%
    filter(spor_start > 0) %>%
    ggplot(aes(x = Year, y = spor_dur)) +
    theme_bw() +
    geom_boxplot(width= .3,alpha = .3, color = "gray") +
    geom_text(aes(year_x_lab, -2),
              label=paste0("italic(p)==",scales::pvalue(as.numeric(anova(m.dur)[[4]][3]))),
              parse=TRUE,
              hjust=1, size=label_size)+
    geom_jitter(height = 0, width = .1)+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank()))
(p23<-
    dtasum %>%
    filter(spor_start > 0) %>%
    ggplot(aes(x = Location, y = spor_dur)) +
    theme_bw() +
    geom_boxplot(width= .3,alpha = .3, color = "gray") +
    geom_text(aes(loc_x_lab, -2),
              label=paste0("italic(p)==",scales::pvalue(as.numeric(anova(m.dur)[[4]][4]))),
              parse=TRUE,
              hjust=1, size=label_size)+
    geom_jitter(height = 0, width = .1)+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank()))

(p31<-
    dtasum %>%
    ggplot(aes(x = Management, y = spor)) +
    theme_bw() +
    geom_boxplot(width= .3,alpha = .3, color = "gray") +
    geom_jitter(height = 0, width = .1)+
    geom_text(aes(2.5, -1),
              label=paste0("italic(p)",
                           scales::pvalue(as.numeric(drop1(m.spor, test = "Chisq")[[5]][3]))),
              parse=TRUE,
              hjust=1, size=label_size
    )+
    ylab("Sporulation intensity(ausc)"))

(p32<-
    dtasum %>%
    ggplot(aes(x = Year, y = spor)) +
    theme_bw() +
    geom_boxplot(width= .3,alpha = .3, color = "gray") +
    geom_jitter(height = 0, width = .1)+
    geom_text(aes(year_x_lab, -1),
              label=paste0("italic(p)==",
                           scales::pvalue(as.numeric(drop1(m.spor, test = "Chisq")[[5]][2]))),
              parse=TRUE,
              hjust=1, size=label_size
    )+
    theme(      axis.title.y=element_blank(),
                axis.text.y=element_blank(),
                axis.ticks.y=element_blank()))

(p33<-
    dtasum %>%
    ggplot(aes(x = Location, y = spor)) +
    theme_bw() +
    geom_boxplot(width= .3,alpha = .3, color = "gray") +
    geom_jitter(height = 0, width = .1)+
    geom_text(aes(loc_x_lab, -1),
              label=paste0("italic(p)",
                           scales::pvalue(as.numeric(drop1(m.spor, test = "Chisq")[[5]][4]))),
              parse=TRUE,
              hjust=1, size=label_size
    )+
    theme(
      axis.title.y=element_blank(),
      axis.text.y=element_blank(),
      axis.ticks.y=element_blank()))




ggpubr::ggarrange(
  p1,p2,p3, 
  p21,p22,p23, 
  p31,p32,p33,
  heights = c(1,1,1.3),
  widths = c(.9,1,1),
  ncol = 3, 
  nrow = 3)
ggsave(here("out/Spor_start_Dur_Intensity.png"), width = 8.3, height = 7)
shell.exec(here("out/Spor_start_Dur_Intensity.png"))


