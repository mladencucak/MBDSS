

list.of.packages <-
  c("tidyverse",
    "readxl",
    "maps",
    "here",
    "ggthemes",
    "ggrepel",
    "ggspatial",
    "sf")
new.packages <-
  list.of.packages[!(list.of.packages %in% installed.packages()[, "Package"])]
#Download packages that are not already present in the library
if (length(new.packages))
  install.packages(new.packages)
if (length(new.packages))
  install.packages(new.packages, repos = c(CRAN = "https://cran.r-project.org/"))

packages_load <-
  lapply(list.of.packages, require, character.only = TRUE)

# if(!"albersusa"%in%installed.packages())  devtools::install_github("hrbrmstr/albersusa")

#Print warning if there is a problem with installing/loading some of packages
if (any(as.numeric(packages_load) == 0)) {
  warning(paste("Package/s", paste(list.of.packages[packages_load != TRUE]), "not loaded!"))
} else {
  print("All packages were successfully loaded.")
}
rm(packages_load, list.of.packages, new.packages)


################################################
#Mapping
################################################
#
#
# browseURL(as.character("https://cran.r-project.org/web/packages/sf/vignettes/sf5.html"))  
# browseURL(as.character("https://www.google.rs/maps/place/WSU+Mt+Vernon+Research+%26+Ext/@45.3951718,-124.0097808,6.25z/data=!4m5!3m4!1s0x54856e5e36e7df0d:0xedebf02c10cde2fa!8m2!3d48.4396286!4d-122.3865128"))  
# 
# #Caclualte distances
# browseURL(as.character("https://ryanpeek.org/2017-10-24-mapping-with-sf-part-1/"))  
# browseURL(as.character("https://gis.stackexchange.com/questions/282750/identify-polygon-containing-point-with-r-sf-package"))  
# 
# #Shape data for states
# browseURL(as.character("https://nceas.github.io/oss-lessons/spatial-data-gis-law/3-mon-intro-gis-in-r.html"))  

df_loc <-
  read.csv(here::here("dat", "disease_loc.csv"), skip = 17) %>% 
  dplyr::select(-starts_with("X"))

df_loc$name <- df_loc$county

df_loc$type <- "inoc"

mtdt <- 
  read_csv( file = here("dat", "wth", "wth_loc_agwet.csv"))


df_loc <-
  mtdt %>% 
  rename(name = stna) %>% 
  mutate(type = "wth") %>% 
bind_rows(df_loc,.)

df_loc_sf <- st_as_sf(
  df_loc,
  coords = c("lon", "lat"),  # for point data
  remove = F, # don't remove these lat/lon cols from df
  crs = 4326) # add projection (this is WGS84)
summary(df_loc_sf)
library(maps)

usa = st_as_sf(maps::map('state',region = c("WA", "OR"), plot = FALSE, fill = TRUE))

st_crs(df_loc_sf) <-
  st_crs(usa)

#########################################
#Map of weather stations
#########################################

 baseplot <- 
ggplot() +
  geom_sf(data = usa,
          color = "darkgray",
          fill = "white") +
  theme_bw(base_family = "Roboto Condensed",
           base_size = 12) +
  annotation_north_arrow(
    location = "br",
    which_north = "true",
    pad_x = unit(2.8, "in"),
    pad_y = unit(0.25, "in"),
    style = north_arrow_fancy_orienteering
  ) +
  labs(x = "Longitude", y = "Latitude") +
  annotation_scale(location = "bl", 
                   # pad_x = 2,
                   width_hint = 0.4) +
  theme(
    panel.background = element_rect(fill = "lightblue"),
    # strip.background = element_blank(),
    # legend.position = c(.18, .91),
    legend.box.background = element_rect(color = "black", size = .5),
    legend.key = element_rect(colour = "transparent", fill = "white")
  ) 

baseplot+
  geom_sf(data = df_loc_sf[df_loc_sf$name != "corvallis", ],
          aes(
             fill = type,
              shape = type,
               color = type
          ),
          # color = "black",
           # shape = 17,
          size = 1.3) +
  coord_sf(
    xlim = c(-123.4, -121.9),
    ylim = c(48.9759, 47.7),
    expand = FALSE
  ) +
  # scale_color_brewer(palette = "Set1")+
  # scale_fill_brewer(palette = "Set1")+
  geom_text_repel(data = df_loc_sf[df_loc_sf$name != "corvallis", ],
                  aes(x = lon, y = lat, label = name,  group = type ),
                  size = 3.2)+
  scale_color_manual(
    name = "Data:",
    labels = c( "Biological sites","Weather stations"),
    values = c("blue", "black")
  ) +
  scale_fill_manual(
    name = "Data:",
    labels = c( "Biological sites","Weather stations"),
    values = c("blue", "black")
  ) +
  scale_shape_manual(
    name = "Data:",
    labels = c( "Biological sites","Weather stations"),
    values = c(16, 17)
  )+
  theme(
    text = element_text(size = 11, family = "TT Times New Roman"),
    legend.position = c(.20, .89) #place legend inside the plotting area
  )

ggsave(
    file = here::here("out" , "map_wth.png"),
    width = 12,
    height = 16,
    units = "cm",
    dpi = 400
  )

shell.exec(here::here("out" , "map_wth.png"))



