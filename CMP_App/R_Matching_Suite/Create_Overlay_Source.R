# Function to create an overlay of two images
Create_Overlay <- function(MapR, Fundus, Match, Replace = FALSE, LimMap = c(200, 400), Pedestal = TRUE) {
  
  # Map the MapR image to the Fundus image using the SLO_to_CMP_Maps function
  MapR <- SLO_to_CMP_Maps(MapR, Fundus, Match)
  
  # Overlay the MapR image on the Fundus image using the Overlay_Maps function
  # Replace: Determines whether the MapR image should replace the Fundus image or be overlaid on top of it
  # LimMap: Specifies the limits of the MapR image
  # Pedestal: Determines whether a pedestal should be added to the final image
  Final_Im <- Overlay_Maps(Fundus = Fundus, MapR = MapR, Replace = Replace,
                           LimMap = LimMap, Pedestal = Pedestal)
  
  # Return the final overlaid image
  return(Final_Im)
}