#3D Vowel Space Density Plots with convex hull imaging
#Vargo, Julian (2025). 3D Vowel Space Density Analyzer [R Script].
#Department of Spanish & Portuguese. University of California, Berkeley.

#Path to your dataset:
csv_path <- "C:/Users/julia/Downloads/research/cook_dispersion/finaldata.csv"

# Unit of measurement for density cutoff is StdDev(n/(ln(Hz))^3)
# Std Dev of number of vowels per log-adjusted Hz cubed (i.e. 3D vowel density SD)
density_cutoff <- 1.5

# How fine grain you want to determine your density radii - think of this as pixel size
# Lower numbers are better for very large datasets. Any density analysis on small datasets is not recommended.
radius_setting <- .1

# How many standard deviations of outlier removal you'd like to implement
outlier_threshold <- 3

######################
# OUTLIER REMOVAL AND VOLUMETRIC MEASUREMENT
######################

# install.packages("ggplot2")
# install.packages("tidyverse")
# install.packages("stringr")
# install.packages("geometry")
# install.packages("plotly")
library(ggplot2)
library(tidyverse)
library(stringr)
library(geometry)
library(plotly)

df <- read.csv(csv_path,header=T,stringsAsFactors=TRUE)

#filter csv so the plotting can be done more efficiently
df <- df %>% filter(str_detect(phoneme, "^(A|E|I|O|U)"))
df$phoneme_nostress <- str_remove(df$phoneme, "[0-2]+$")

# log_e normalization of the data before outlier removal
df <- df %>% mutate (n_f1 = log(w_f1))
df <- df %>% mutate (n_f2 = log(w_f2))
df <- df %>% mutate (n_f3 = log(w_f3))

#Outlier removal for vowels
df <- df %>%
  group_by(phoneme) %>%
  mutate(
    f1_centroid = median(n_f1, na.rm = TRUE),
    f2_centroid = median(n_f2, na.rm = TRUE),
    f3_centroid = median(n_f3, na.rm = TRUE),
    distance_from_centroid = sqrt(
      (n_f1 - f1_centroid)^2 +
        (n_f2 - f2_centroid)^2 +    
        (n_f3 - f3_centroid)^2
    ),
    outlier_limit = outlier_threshold * sd(distance_from_centroid, na.rm = TRUE)
  ) %>%
  filter(distance_from_centroid < outlier_limit) %>%
  ungroup()

#Now each vowel must be assigned a density depending on an arbitrarily specified spherical radius.
#The radius for each point will be .1 times the sqrt of std devs of each f1, f2, f3

f1_sd = sd(df$n_f1)
f2_sd = sd(df$n_f2)
f3_sd = sd(df$n_f3)
pythag_sd = sqrt(f1_sd^2 + f2_sd^2 + f3_sd^3)
density_radius = radius_setting * pythag_sd

#Given our radius, we'll calculate the number of points in a sphere
#The way the coding logic works is, we increase a count only if the radius between two points is less than our target radius
#We'll check every single point, so this is computationally intensive
local_density <- numeric(nrow(df))
for (i in 1:nrow(df)) {
  distances_squared <- (df$n_f1 - df$n_f1[i])^2 +
    (df$n_f2 - df$n_f2[i])^2 +
    (df$n_f3 - df$n_f3[i])^2
  local_density[i] <- sum(distances_squared <= density_radius) / ((4/3) * pi * density_radius^3)
}

df$local_density <- local_density
avg_density <- mean(df$local_density)
sd_density <- sd(df$local_density)
df <- df %>% filter(local_density > (avg_density - (1.5*sd_density)))

#Calculate the convex hull of the total dataset
points <- cbind(df$n_f1, df$n_f2, df$n_f3)
points <- as.matrix(points)
hull <- convhulln(points)
hull_volume <- convhulln(points, options = "FA")
paste("The total hull volume is: ", hull_volume$vol)

#convex hull volume for each file name
full_speech <- df %>% filter(df$file_name == "full_speech.TextGrid")
full_speech_points <- cbind(full_speech$n_f1, full_speech$n_f2, full_speech$n_f3)
full_speech_points <- as.matrix(full_speech_points)
full_speech_hull <- convhulln(full_speech_points)
full_speech_hull_volume <- convhulln(full_speech_points, options = "FA")
paste("The full_speech hull volume is: ", full_speech_hull_volume$vol)

homepods <- df %>% filter(df$file_name == "homepods.TextGrid")
homepods_points <- cbind(homepods$n_f1, homepods$n_f2, homepods$n_f3)
homepods_points <- as.matrix(homepods_points)
homepods_hull <- convhulln(homepods_points)
homepods_hull_volume <- convhulln(homepods_points, options = "FA")
paste("The homepods hull volume is: ", homepods_hull_volume$vol)

keynote <- df %>% filter(df$file_name == "keynote.TextGrid")
keynote_points <- cbind(keynote$n_f1, keynote$n_f2, keynote$n_f3)
keynote_points <- as.matrix(keynote_points)
keynote_hull <- convhulln(keynote_points)
keynote_hull_volume <- convhulln(keynote_points, options = "FA")
paste("The keynote hull volume is: ", keynote_hull_volume$vol)

mit_speech <- df %>% filter(df$file_name == "mit_speech.TextGrid")
mit_speech_points <- cbind(mit_speech$n_f1, mit_speech$n_f2, mit_speech$n_f3)
mit_speech_points <- as.matrix(mit_speech_points)
mit_speech_hull <- convhulln(mit_speech_points)
mit_speech_hull_volume <- convhulln(mit_speech_points, options = "FA")
paste("The mit_speech hull volume is: ", mit_speech_hull_volume$vol)

one_on_one <- df %>% filter(df$file_name == "one_on_one.TextGrid")
one_on_one_points <- cbind(one_on_one$n_f1, one_on_one$n_f2, one_on_one$n_f3)
one_on_one_points <- as.matrix(one_on_one_points)
one_on_one_hull <- convhulln(one_on_one_points)
one_on_one_hull_volume <- convhulln(one_on_one_points, options = "FA")
paste("The one_on_one hull volume is: ", one_on_one_hull_volume$vol)

##########################
# VISUALIZATION
##########################
plot_ly(
  data = df,
  x = ~n_f1,
  y = ~n_f2,
  z = ~n_f3,
  type = "scatter3d",
  mode = "markers",
  marker = list(size = 3, color = "blue", opacity = .2)
)

x <- df$n_f1
y <- df$n_f2
z <- df$n_f3
i <- hull[,1] - 1
j <- hull[,2] - 1
k <- hull[,3] - 1

plot_ly() %>%
  add_markers(
    data = df,
    x = ~n_f1,
    y = ~n_f2,
    z = ~n_f3,
    marker = list(size = 2, color = "blue", opacity = 0.3)
  ) %>%
  add_mesh(
    x = x,
    y = y,
    z = z,
    i = i,
    j = j,
    k = k,
    intensity = rep(1, length(x)),
    colorscale = list(c(0, "lightgrey"), c(1, "lightgrey")),
    opacity = 0.3,
    showscale = FALSE
  )