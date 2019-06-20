
############################# pERPred hex sticker #############################

library(here)
library(tidyverse)
library(hexSticker)

sticker(
  here("transparent_brain.png"),
  package = "pERPred", # package name
  p_size = 8,          # font size for package name
  p_x = 1,             # x position for package name (1 = center)
  p_y = 0.7,           # y position for package name (1 = center)
  s_x = 1,             # x position for plot (1 = center)
  s_y = 1.15,          # y position for plot (1 = center)
  asp = 0.8,           # aspect ratio for plot
  s_width = 0.55,       # width for plot
  h_fill = "#22AA9A",  # main color
  h_color = "#B8C9C7", # border color
  l_x = 1,
  l_y = 1,
  l_alpha = 0.2,
  l_height = 4,
  l_width = 4,
  spotlight = TRUE
)
