## Discrete Origin Assignment of Wisconsin Mallards - WI/MN DNR Isotope Project
## Penelope Murphy, code adapted from Drew Fowler

rm(list=ls()) 

################################
## PACKAGES
if (!require('raster')) install.packages('raster'); library('raster')
if (!require('terra')) install.packages('terra'); library('terra')
if (!require('sf')) install.packages('sf'); library('sf')
if (!require('sp')) install.packages('sp'); library('sp')
if (!require('rnaturalearth')) install.packages('rnaturalearth'); library('rnaturalearth')
if (!require('rnaturalearthdata')) install.packages('rnaturalearthdata'); library('rnaturalearthdata')
if (!require('rnaturalearthhires')) install.packages('rnaturalearthhires'); library('rnaturalearthhires')
if (!require('dplyr')) install.packages('dplyr'); library('dplyr')
if (!require('readr')) install.packages('readr'); library('readr')
if (!require('RColorBrewer')) install.packages('RColorBrewer'); library('RColorBrewer')
if (!require('ggplot2')) install.packages('ggplot2'); library('ggplot2')
if (!require('viridis')) install.packages('viridis'); library('viridis')
if (!require('rmapshaper')) install.packages('rmapshaper'); library('rmapshaper')
if (!require('tidyr')) install.packages('tidyr'); library('tidyr')
if (!require('lubridate')) install.packages('lubridate'); library('lubridate')

##############################

## PRECIPITATION ISOSCAPE

# download precip isoscape from github repository
download.file(
  "https://raw.githubusercontent.com/pwmurphy/Isotope_Discrete_Assignment_Demo/main/data/d2h_GS.tif",
  destfile = "d2h_GS.tif",
  mode = "wb"
)

# load precipitation isoscape raster and set CRS
isoscape_GS <- rast("d2h_GS.tif")
crs(isoscape_GS) <- "+proj=longlat +proj=longlat +datum=WGS84 +no_defs"

plot(isoscape_GS) #view precipitation isoscape

## FEATHER ISOSCAPE CALIBRATION

# create feather isoscape using the mean growing season (MGS) dabbler calibration equation from Kusack et al. 2023
isoscape_GS_NA <- (-69.9) + (0.7 * isoscape_GS)

plot(isoscape_GS_NA) # view isoscape


## DEFINE ORIGIN REGION VALUE BREAKPOINTS
# To define origin region value breakpoints, use the minimum and maximum isoscape values for Wisconsin

# load in US states
state_basemap <- ne_states(country = "united states of america", # load in all the US states
                           returnclass = "sf") 

# isolate Wisconsin
wi_basemap <- state_basemap %>%
  filter(fips == 'US55') %>% # select Wisconsin
  st_transform("+proj=longlat +proj=longlat +datum=WGS84 +no_defs") %>% # set crs
  ms_simplify(keep=0.05) %>% # simplify the polygon shape a bit
  as("Spatial") %>% 
  vect() # makes it a SpatVector to work with the raster

# mask and crop the feather isoscape to Wisconsin
wi_isoscape_NA <- mask(isoscape_GS_NA, wi_basemap)
wi_isoscape_NA <- crop(wi_isoscape_NA,wi_basemap)

# plot isoscape
plot(wi_isoscape_NA)

# nice ggplot of Wisconsin isoscape map, for saving if desired
# plot raster
wi_isoscape_NA_df <- as.data.frame(wi_isoscape_NA, xy = TRUE, na.rm = TRUE) %>% 
  rename(value = d2h_GS)

ggplot() +  
  geom_tile(data=wi_isoscape_NA_df, aes(x=x, y=y, fill=value), alpha=0.8) + 
  scale_fill_viridis(name = expression(paste(delta^{2}, "H (\u2030)"))) +
  coord_equal() +
  theme_minimal()+
  theme(legend.position="right") +
  theme(legend.key.width=unit(2,units = "cm" ))+
  labs(
    x = NULL,
    y = NULL,
    title = 'Wisconsin Feather Isoscape')

# Isolate range of WI feather d2h values
wi_d2h_values_NA <- as.data.frame(values(wi_isoscape_NA)) %>%
  drop_na() %>% 
  rename(d2h_vs_vsmow_original = d2h_GS) %>% 
  mutate(source = "wi_GS_feather_isoscape")

# establish feather d2h value boundaries for local versus non local origin using minimum and maximum expected d2h values in WI +/- the standard deviation of the calibration equation
SD_eq <- 17.7 # taken from Kusack et al. 2023 MGS Dabbler Calibration Equation

local_max <- max(wi_d2h_values_NA$d2h_vs_vsmow_original) + SD_eq
local_min <- min(wi_d2h_values_NA$d2h_vs_vsmow_original) - SD_eq

## DISCRETE ASSIGNMENT OF 'LOCAL' VS 'NORTHERN' OR 'SOUTHERN' ORIGIN

# download feather d2h sample results for after hatch-year Mallards harvested in Wisconsin, unless you cloned the repo to your local machine
download.file(
  "https://raw.githubusercontent.com/pwmurphy/Isotope_Discrete_Assignment_Demo/main/data/MALL_WI_ahy_dataset.csv",
  destfile = "MALL_WI_ahy_dataset.csv",
  mode = "wb"
)

# load in feather dataset
mall_feather_df <- read_csv("MALL_WI_ahy_dataset.csv") %>% 
  mutate(collection_date = mdy(collection_date)) # specify collection_date as a date

# use the range of expected d2h feather values in Wisconsin to assign local or northern origin "predicted origin" to each sample, amend to dataframe 
mall_feather_df <- mall_feather_df %>% 
  mutate(predicted_origin = case_when(d2h_vs_vsmow_original <= local_min ~ "northern", # more negative values are 
                                      d2h_vs_vsmow_original < local_max & d2h_vs_vsmow_original > local_min ~ "local",
                                      d2h_vs_vsmow_original >= local_max ~ "southern"),
         doy = yday(collection_date))


# Figure of results: number of birds assigned as 'local', 'nothern' or 'southern' by discrete assignment
mall_ahy_sample_size <- count(mall_feather_df %>% drop_na(d2h_vs_vsmow_original)) # pull sample size

mall_feather_df %>% drop_na(d2h_vs_vsmow_original) %>%
  group_by(sex, predicted_origin) %>%
  tally() %>%
  mutate(freq = n / sum(n)) %>%
  ggplot(aes(x = sex, y = n, fill = predicted_origin )) +
  geom_col(position = "dodge") +
  theme_minimal()+
  geom_text(position = position_dodge(width = 0.9),
            hjust = -0.25, aes(label = scales::percent(freq, accuracy = 0.1))) +
  scale_fill_manual(values = c("local" = "black","northern" = "#9fb052", "southern" = "#9988bd"))+
  ylim(c(0, 100)) +
  coord_flip() +
  labs(
    y = 'Number of individuals',
    x = '',
    title = paste0('After Hatch-year Wisconsin-harvested Mallards: \nDiscrete isotope assignment (n = ', mall_ahy_sample_size,')'),
    fill = "Predicted origin") 


## CONSIDERATION: IS DISCRETE ASSIGNMENT TOO COARSE?

# To determine what area corresponds to our 'local' origin classification, 
# we can look at the feather isoscape across the breeding range of the Mississippi and Central Flyways 
# to see what the 'local' origin range covers. 

# Download mallard breeding (MS and CN flyway) assignment range from Github repo location
download.file(
  "https://raw.githubusercontent.com/pwmurphy/Isotope_Discrete_Assignment_Demo/main/data/mall.assignment.range.rds",
  destfile = "mall.assignment.range.rds",
  mode = "wb"
)

# crop and mask isoscape raster to assignment range and set CRS
mall.assignment.range <- readRDS(file = "mall.assignment.range.rds") %>% 
  vect(crs = isoscape_GS_NA)

range_isoscape <- crop(isoscape_GS_NA, mall.assignment.range)
range_isoscape <- mask(range_isoscape, mall.assignment.range) # warning about CRS mismatch can be ignored: a difference in 'remarks'

plot(range_isoscape)

# upload country, state/province shapefiles

northamerica <- ne_countries(continent = "North America", scale = 50, returnclass = "sf") %>% # Countries
  st_transform(st_crs('EPSG:4326')) # Project to WGS84

northamerica.states <- ne_states(country =  c("Canada","United States of America"), returnclass = "sf") %>% # States and provinces
  st_transform(st_crs('EPSG:4326')) #%>% # Project to WGS84

plot(st_geometry(northamerica.states), add = T)

# feclassify the raster surface according to the discrete origin region breakpoints to determine which part of the flyways 
# are considered 'northern', 'local', and 'southern'

## create a matrix of the discrete origin area value ranges
reclass_matrix_dabblers <- matrix(c(-Inf, local_min, 0, # northern origin area
                                    local_min, local_max, 1, # local origin area
                                    local_max, Inf, 2), # southern origin area
                                  ncol = 3, byrow = TRUE)


# Reclassify the raster with the matrix
range_isoscape_binary <- classify(range_isoscape, reclass_matrix_dabblers)
category_colors <- c("#9fb052", "#9988bd", "black")

# make nice plot of the flyway feather isoscape dabblers, 'local' area

# Prepare the reclassified isoscape for plotting in ggplot
range_binary_df <- as.data.frame(range_isoscape_binary, xy = TRUE, na.rm = FALSE) %>% 
  rename(value = d2h_GS) %>% 
  mutate(region = case_when(
    value == 0 ~ "northern",
    value == 1 ~ "local",
    value == 2 ~ "southern",
    is.na(value) ~ "out of range"
  ))

# Get the extent of the raster to use for the map extent
range_extent <- ext(range_isoscape_binary)

# Create the plot
ggplot() +
  # Plot the raster
  geom_tile(data = range_binary_df, aes(x = x, y = y, fill = as.factor(region))) +
  scale_fill_manual(values = c("local" = "black","northern" = "#9fb052", "southern" = "#9988bd", "out of range" = "white")) +
  # Plot the sf dataframe
  geom_sf(data = northamerica.states, fill = NA, color = "black", linewidth = 0.5) +
  # Set the coordinate limits to match the raster extent
  coord_sf(
    xlim = c(range_extent$xmin, range_extent$xmax),
    ylim = c(range_extent$ymin, range_extent$ymax),
    expand = FALSE
  ) +
  # Add labels and theme
  labs(
    title = "Discrete assignment origin regions for Dabblers",
    x = "Longitude",
    y = "Latitude",
    fill = "Discrete origin region"
  ) +
  theme_minimal()


##  EXTRA: ASSESS HOW TIMING OF HARVEST IS RELATED TO ORIGIN LOCATION 

# Plot after hatch-year mallard d2H by feather value
mall_feather_df %>% 
  drop_na(d2h_vs_vsmow_original) %>%
  ggplot(aes(x=doy, y=d2h_vs_vsmow_original)) +
  geom_point() +
  geom_smooth(method=lm , fill="grey", se=T) +
  theme_minimal()+
  ylab(expression(paste(delta^{2}, "H (\u2030)")))+
  xlab("Julian day")+
  ggtitle(expression(paste("After hatch-year WI harvest mallard feather ", delta^{2}, "H")))


# Scatter plot of AHY mallard d2H by latitude and by sex and predicted origin 
mall_feather_df %>% drop_na(d2h_vs_vsmow_original) %>%
  ggplot(aes(x=doy, y=d2h_vs_vsmow_original, color=predicted_origin, shape=sex )) +
  geom_point() +
  geom_smooth(method=lm , fill="gray", se=T) +
  theme_minimal() +
  scale_color_manual(values = c("local" = "black","northern" = "#9fb052", "southern" = "#9988bd")) +
  labs(
    color = "Predicted origin",
    shape = "Sex",
    y = expression(paste(delta^{2}, "H (\u2030)")),
    x = "Julian day",
    title = expression(paste("After hatch-year WI harvest mallard feather ", delta^{2}, "H)"))
  )


sessionInfo()
