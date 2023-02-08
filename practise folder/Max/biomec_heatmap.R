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
cols <- tail(p_data$X, n = 46)

        
plot <- ggplot(data = p_data, aes(X, Y, fill = Z))+
        geom_tile()+
        scale_color_viridis_d(n)+
        transition_time(frames)+
        ease_aes('linear')

setwd("/Users/maxul/Documents/Skole/Master 21-22/From Naomi/")


animate(plot)

anim_save("pressure_plateC.gif", animation = last_animation(), path = "/Users/maxul/Documents/Skole/Master 21-22/From Naomi/")



### make dataframe of x and y coordinate for every newton of force

p_data %>% 
        group_by(frames) %>% 
        summarise(Z = sum(Z)) %>% 
        filter(Z == max(Z))

p_data %>% 
        filter(frames == "17") %>% 
        select(1:3) -> max_frame



max_frame[1,1:2]

x_y <- data.frame(matrix(NA, nrow = 1, ncol = 2,))

colnames(x_y) <- c("Y", "X")

for (i in 1:nrow(max_frame)) {
 
        if (max_frame[i,3] > 0) {
                x_y[nrow(x_y)+1,] <- as.data.frame(lapply(max_frame[i,1:2], rep, times =  max_frame[i,3]))
                 
                
        }
     
      print(i)
}

x_y <- na.omit(x_y)

# remove "Col" from Y value

x_y[,2] <- gsub("Col", "", x_y$Y)


### run kmeans clustering on x_y

set.seed(6)
X <- x_y[1:2]

xyss <- vector()

for (i in 1:10) {
        xyss[i] <- sum(kmeans(X,i)$withinss)
}
plot(1:10, xyss, type = "b", main = "clusters of force on plate", xlab = "number of clusters", ylab = "XYSS")


# make cluster with identified elbow
set.seed(29)
kmeans <- kmeans(X, 10, iter.max = 300, nstart = 10)



# visualize clusters

library(cluster)

clusplot(X, kmeans$cluster, lines = 0, shade = TRUE, color = TRUE, labels = 2, plotchar = FALSE, span = TRUE, main = "clusters of force on plate", xlab = "number of clusters", ylab = "XYSS" )


# add cluster number to data frame

xy_cluster <- x_y %>% 
        mutate(cluster = kmeans$cluster)

xy_cluster %>% 
        ggplot(aes(x = X, y = Y, color = as.factor(cluster)))+
        geom_count()

################################

### repeat with all datapoints


x_y2 <- data.frame(matrix(NA, nrow = 1, ncol = 2,))

colnames(x_y2) <- c("X", "Y")

p_data2 <- p_data %>% 
        filter(Z >0)

for (i in 1:nrow(p_data2)) {
        
        if (p_data2[i,3] > 0) {
                x_y2[nrow(x_y2)+1,] <- as.data.frame(lapply(p_data2[i,1:2], rep, times =  p_data2[i,3]))
                
                
        }
        
        print(i)
}






for (i in 1:nrow(p_data2)) {
        
        for (j in 1:(p_data2$Z)) {
                x_2[nrow(x_y2)+1,] <- as.data.frame(p_data2[i, 1:2])
        }
                
                
        
        print(i)
}


x_y <- data.frame()

for (i in 1:nrow(p_data2)) {
        repeat_times <- p_data2[i,3]
        repeat_rows <- data.frame(p_data2[i,1:2])
        repeat_rows[nrow(repeat_rows)+1,] <- as.data.frame(lapply(p_data2[i,1:2], rep, times =  repeat_times))
        x_y <- rbind(x_y, repeat_rows)
}

x_y<- na.omit(x_y)

# make cluster with identified elbow
set.seed(29)
x2 <- x_y %>% 
        mutate(X = gsub("Col", "", X))


kmeans2 <- kmeans(x = x2, centers = 5, iter.max = 300, nstart = 10)


xy_cluster <- x2 %>% 
        mutate(cluster = kmeans2$cluster)

xy_cluster %>% 
        ggplot(aes(x = X, y = Y, color = as.factor(cluster)))+
        geom_count()+
        scale_size_binned_area(breaks = scales::extended_breaks())

xy_cluster



