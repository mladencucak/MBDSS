



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
    "conflicted",
    "lubridate",
    "naniar",
    "tibble",
    "tidyr",
    "glmmTMB",
    "plotly"
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
conflict_prefer("layout", "plotly")

rm(packages_load, list.of.packages, new.packages)


if(!c("insight") %in% installed.packages())install.packages(c("insight"), repos = 'https://easystats.r-universe.dev')
if(!c("insight") %in% installed.packages())install.packages(c("performance"), repos = 'https://easystats.r-universe.dev')

                                                            
lapply(c("insight", "performance"), require, character.only = TRUE)


 #first we register to upload to plotly
# Please register 
Sys.setenv("plotly_username"="XXXX")
Sys.setenv("plotly_api_key"="XXXX")


##################################################################
#Original data
##################################################################
dis_df <- read_delim(here("scr", "model", "par", "inf_period_hildebrand.txt"), delim = ",")


dis_df$wet_dur <- round(dis_df$wet_dur,0)
dis_df$dis_prop2 <- dis_df$dis_prop/100
dis_df$dis_prop2[dis_df$dis_prop2 == 1] <- .999
dis_df$dis_prop2[dis_df$dis_prop2 == 0] <- .001

tit <- "Hildebrandt figure"
ggplot(dis_df,aes(x= wet_dur, y = dis_prop, group = factor(temp), colour = factor(temp)))+
  geom_point()+
  geom_line()+
  scale_color_brewer(palette="Dark2")+
  labs(
    title = tit,
    x = "Wetness duration (h)",
    y = "Infection (%)",
    color = "Temperature (°C)",
    caption = expression(paste(
      "\n  "
    )))+
  egg::theme_article()+
  scale_y_continuous(limits = c(0,100))+
  theme(legend.position = "top")+
  theme(
    axis.text.x = element_text(size = 10),
    # axis text size
    axis.text.y = element_text(vjust = 0.2),
    # axis text alignment
    axis.ticks = element_line(size = 0.4),
    axis.title = element_text(size = 11, face = "bold"),
    # axis title size and bold
  ) 

ggsave(
  filename = here::here("scr", "model", "fig", paste0(tit, ".png")),
  width = 8.2,
  height = 7,
  dpi = 300
)


ggplot(dis_df,aes(x= temp, y = dis_prop, group = factor(wet_dur), colour = factor(wet_dur)))+
  geom_point()+
  geom_line()+
  scale_color_brewer(palette="Dark2") 



##################################################################
# Fit model
##################################################################
 

dis_df <- 
  bind_rows(
    dis_df,
    dis_df[dis_df$temp == 2, ]%>% 
      mutate( temp =30)
  )


 

fit_noint <- glmmTMB(dis_prop2 ~ poly(temp,3) + (wet_dur + I(log(wet_dur+1))), 
                       family = beta_family, 
                       data = dis_df)
fit_int   <- glmmTMB(dis_prop2 ~ poly(temp,3) * (wet_dur + I(log(wet_dur+1))), 
                             family = beta_family, 
                             data = dis_df)
anova( fit_noint, fit_int)



summary(fit_noint)
summary(fit_int)
glmmTMB:::Anova.glmmTMB(fit_noint)
glmmTMB:::Anova.glmmTMB(fit_int)

performance::r2(fit_noint)
performance::r2(fit_int)


fit <- fit_int

df_fit <- cbind(dis_df,fit =plogis(predict(fit,dis_df))*100) 
df_fit
str(fit)
conf_int <- predict(fit, se.fit = TRUE)
upr <- plogis(conf_int$fit + 1.96 * conf_int$se.fit)*100
lwr <- plogis(conf_int$fit - 1.96 * conf_int$se.fit)*100

ggplot(data = df_fit,aes(x= wet_dur, y = dis_prop, colour = "Observed"))+
  geom_point()+
  geom_line()+
  geom_line(aes(wet_dur, fit, colour = "Predicted")) +
  geom_ribbon(aes(ymin=lwr,ymax=upr), alpha=0.3)+
  scale_color_manual(values = c( 
    Observed = "black", 
    Predicted = "red")) + 
  ylim(c(0,100))+
  xlab("Wetness Duration (h)")+
  ylab("Disease severity (%)")+
  facet_grid(~temp)+
  theme_bw()

ggsave(
  filename = here::here("scr", "model", "fig", paste0("wet_dur_fit_",deg, ".png")),
  width = 8.2,
  height = 7,
  dpi = 300
)

ggplot(data = df_fit,aes(x= temp, y = dis_prop, colour = "Observed"))+
  geom_point()+
  geom_line()+
  geom_line(aes(temp, fit, colour = "Predicted")) +
  geom_ribbon(aes(ymin=lwr,ymax=upr), alpha=0.3)+
  scale_color_manual(values = c( Observed = "black", Predicted = "red"))+
  facet_grid(~wet_dur, scales = "free")+
  xlab("Temperature (˚C)")+
  ylab("Disease severity (%)")+
  theme_bw()
ggsave(
  filename = here::here("scr", "model", "fig", paste0("temp_fit_",deg, ".png")),
  width = 8.2,
  height = 7,
  dpi = 300
)
shell.exec(here::here("scr", "model", "fig"))

mod <- fit
save(mod, file = here("scr/model/inf_model.RData"))

x.seq <- seq(min(dis_df$temp, na.rm=TRUE), max(dis_df$temp, na.rm=TRUE), by = 1) 
y.seq <- seq(min(dis_df$wet_dur, na.rm=TRUE), max(dis_df$wet_dur, na.rm=TRUE), by = 1)

predfun <- function(x,y){
  newdat <- data.frame(temp = x, wet_dur=y)
  plogis(predict(fit, newdata=newdat))*100
}

fit <- outer(x.seq, y.seq, Vectorize(predfun))

fit[ ,1:2] <- 0
 

f1 <- list(
  family = "Arial, sans-serif",
  size = 16,
  color = "grey"
)
f2 <- list(
  family = "Old Standard TT, serif",
  size = 12,
  color = "black"
)
tit<- list(
  title = "Infection Function",
  titlefont = f1,
  tickfont = f2
)
x <- list(
  title = "Temperature (˚C)",
  titlefont = f1,
  showticklabels = TRUE,
  tickangle = 30,
  tickfont = f2,
  exponentformat = "E",
  nticks= length(unique(dis_df$temp)),
  range= c(min(x.seq), max(x.seq))%>% rev(),
  tickvals =unique(dis_df$temp),
  ticktext = unique(dis_df$temp),
  backgroundcolor="white",
  gridcolor="black",
  showbackground=TRUE,
  zerolinecolor="black"
)
y <- list(
  title = "Wetness Duration (h)",
  titlefont = f1,
  showticklabels = TRUE,
  tickangle = -45,
  tickfont = f2,
  exponentformat = "E",
  nticks = length(unique(dis_df$wet_dur)),
  range= c(min(y.seq), max(y.seq))%>% rev(),
  tickvals =unique(dis_df$wet_dur),
  ticktext = unique(dis_df$wet_dur),
  backgroundcolor="white",
  gridcolor="black",
  showbackground=TRUE,
  zerolinecolor="black"
)
z <- list(
  title = "Infection risk (%)",
  titlefont = f1,
  showticklabels = TRUE,
  tickvals =seq(0.1,1, .1),
  ticktext = paste(seq(10,100, 10)),
  tickangle = 0,
  tickfont = f2,
  range = c(0,1),
  # exponentf
  zerolinecolor="black"
)



# ##Custom ticks
# axx <- list(
#   ticketmode = 'array',
#   ticktext = c("Huey", "Dewey", "Loormat = "E",
#   nticks= 10,
#   range= c(0, 100),
#   ticktext = paste(seq(0,10, 10)),
#    backgroundcolor="white",
#   gridcolor="black",
#   showbackground=TRUE,uie"),
#   tickvals = c(0,25,50),
#   range = c(-25,75)
# )


(p <-
plot_ly() %>%
  add_markers(
    x = ~ dis_df$temp,
    y = dis_df$wet_dur,
    z = dis_df$dis_prop/100,
    symbol = 'triangle',
    marker = list(
      color = ~ dis_df$dis_prop/100,
      size = 5,
      line = list(color = "black",
                  width = 1)
    )
  ) %>%
  add_surface(
    x = ~ x.seq,
    y = ~ y.seq,
    z = t(fit/100),
    opacity = .6
  ) %>%
  layout(title = tit,
         scene = list(
           xaxis = x,
           yaxis = y,
           zaxis = z
         )))



plotly_IMAGE(p, width = 750, height = 550, format = "png", scale = 2,
             out_file = here("scr/model/fig/surface_inf.png"))
shell.exec(here("scr/model/fig/surface_inf.png"))




