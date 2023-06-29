moving_average = function(y, window_length){
  T_data = length(y)
  y_output = y
  if (window_length %% 2 != 0){  # odd window_length
    temp1 = round(window_length/2)
    temp2 = T_data - round(window_length/2) + 1
    half_length = (round(window_length/2)-1)
    for (i in temp1:temp2){
      window_start =  i - half_length
      window_end = i + half_length
      y_output[i] = mean(y[window_start:window_end])
      
    }
  }
  
  return(y_output)
}

plot_graphs = function(x, adj_block){
  # plot graph on a world map
  # x: locations
  # graph
  p = dim(adj_block)[1]
  require(grid)
  center <- 180
  # shift coordinates to recenter x
  long.recenter <- ifelse(x[,1] < center - 180 , x[,1] + 360, x[,1])
  x = data.frame(long.recenter, latitude = x[,2])
  # shift coordinates to recenter worldmap
  worldmap <- map_data ("world", wrap = c(0, 360))
  # Plot worldmap using data from worldmap.cp
  plt = ggplot(aes(x = long, y = lat), data = worldmap)+
    coord_cartesian(xlim = c(180-80, 180+110), ylim = c(-50, 70))+
    geom_polygon(aes(group = group), fill="#f9f9f9", colour = "grey65") + 
    theme_bw() +
    geom_point(data = x,
               aes(x = long.recenter, y = latitude),
               pch = 19, size = 2, alpha = .4, col = 'blue') +
    xlab("") + ylab("")
  
  # bottom left to top right | top right to bottom left
  myCurve1<-curveGrob(0, 0, 1, 1, default.units = "npc",
                      curvature = -0.3, angle = 60, ncp = 40, shape = 1,
                      square = FALSE, squareShape = 1,
                      inflect = FALSE, open = TRUE,
                      debug = FALSE,
                      name = NULL, gp = gpar(), vp = NULL)
  # botton right to top left | top left to bottom right
  myCurve2<-curveGrob(0, 1, 1, 0, default.units = "npc",
                      curvature = -0.3, angle = 60, ncp = 40, shape = 1,
                      square = FALSE, squareShape = 1,
                      inflect = FALSE, open = TRUE,
                      debug = FALSE,
                      name = NULL, gp = gpar(), vp = NULL)
  # horizontal
  myCurve3<-curveGrob(0, 0, 1, 0, default.units = "npc",
                      curvature = -0.3, angle = 60, ncp = 40, shape = 1,
                      square = FALSE, squareShape = 1,
                      inflect = FALSE, open = TRUE,
                      debug = FALSE,
                      name = NULL, gp = gpar(), vp = NULL)
  # verticle
  myCurve4<-curveGrob(0, 1, 0, 0, default.units = "npc",
                      curvature = -0.3, angle = 60, ncp = 40, shape = 1,
                      square = FALSE, squareShape = 1,
                      inflect = FALSE, open = TRUE,
                      debug = FALSE,
                      name = NULL, gp = gpar(), vp = NULL)
  
  
  for (j_1 in 1:(p-1)){
    for (j_2 in (j_1+1):p){
      if (adj_block[j_1,j_2]){
        temp1 = ((x[j_1,1] < x[j_2,1]) & (x[j_1,2] < x[j_2,2])) | 
          ((x[j_1,1] > x[j_2,1]) & (x[j_1,2] > x[j_2,2]))
        temp2 = ((x[j_1,1] < x[j_2,1]) & (x[j_1,2] > x[j_2,2])) | 
          ((x[j_1,1] > x[j_2,1]) & (x[j_1,2] < x[j_2,2]))
        temp3 = x[j_1,2] == x[j_2,2] 
        temp4 = x[j_1,1] == x[j_2,1]
        if (temp1){
          plt <- plt +  
            annotation_custom(grob=myCurve1,x[j_1,1],x[j_2,1],x[j_1,2],x[j_2,2])
        }else if (temp2){
          plt <- plt +  
            annotation_custom(grob=myCurve2,x[j_1,1],x[j_2,1],x[j_1,2],x[j_2,2])
        }else if ( temp3 ){
          plt <- plt +  
            annotation_custom(grob=myCurve3,x[j_1,1],x[j_2,1],x[j_1,2],x[j_2,2])
        }else if ( temp4 ){
          plt <- plt +  
            annotation_custom(grob=myCurve4,x[j_1,1],x[j_2,1],x[j_1,2],x[j_2,2])
        }else{
          print(c(x[j_1,1],x[j_2,1],x[j_1,2],x[j_2,2]))
        }
      }
    }
  }
  return(plt)
  
}
