#----------------------------------------------------------------
#
#   plotting pressure data heatmap
#
#
#----------------------------------------------------------------


library(ggplot2); library(tidyverse); library(readxl); library(gganimate); library(gifski)


sheets <- excel_sheets("C:/Users/maxul/Documents/Skole/Master 21-22/From Naomi/FramePressureData.xlsx")

pressure_df <- list()

for (i in 1:length(sheets)) {
        x <- read_excel("C:/Users/maxul/Documents/Skole/Master 21-22/From Naomi/FramePressureData.xlsx", sheet = i)
        pressure_df[[i]] <- as.data.frame(x)
        print(i)
}

pressure_data = list()

for( i in 1: length(pressure_df)){
        pressure_df[i] %>% 
                as.data.frame() %>% 
                mutate(Y = rownames(as.data.frame(pressure_df[i]))) %>% 
                pivot_longer(names_to = "X", values_to = "Z", cols = 1:96) %>% 
                mutate(frames = (i))-> t
        pressure_data[[i]]<- t
        
        
}
cols <- tail(p_data$X, n = 46)

p_data <- bind_rows(pressure_data) %>% 
        filter(Y %in% 20:66,
               X %in% cols)


        
plot <- ggplot(data = p_data, aes(X, Y, fill = Z))+
        geom_tile()+
        scale_color_viridis_d(n)+
        transition_time(frames)+
        ease_aes('linear')

setwd("/Users/maxul/Documents/Skole/Master 21-22/From Naomi/")


animate(plot)

anim_save("pressure_plateC.gif", animation = last_animation(), path = "/Users/maxul/Documents/Skole/Master 21-22/From Naomi/")
