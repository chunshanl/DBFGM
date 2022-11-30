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