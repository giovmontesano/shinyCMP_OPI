#Checks if polygon is clockwise (necessary for later)
clockwise <- function(x, y) {
  
  x.coords <- c(x, x[1])
  y.coords <- c(y, y[1])
  
  double.area <- sum(sapply(2:length(x.coords), function(i) {
    (x.coords[i] - x.coords[i-1])*(y.coords[i] + y.coords[i-1])
  }))
  
  double.area > 0
} 


## Snell law ######
Snell_law <- function(Normal_Axis, mInc, n1, n2){
  
  Interface_Angle <- pi/2 + Normal_Axis;
  mN <- tan(Normal_Axis);
  mI <- tan(Interface_Angle);
  
  IncidentAngle <- ((atan((mN - mInc)/(1 + mN*mInc))));
  RefAngle <- -asin(sin(IncidentAngle)*(n1/n2));
  
  mRef <- tan(RefAngle + Normal_Axis);
}

## Schematic eye ########
schem_eye <- function(Axial_Length = (11.459/11.06*23.01), ExpMod = 'Global', K = 43.0769, alpha = 0){
  
  Schematic_Eye <- data.frame(C_Index = 1.336, A_Index = 1.336, 
                              V_Index = 1.336, Lens_Index = 1.43, Axial_Length = Axial_Length, ExpMod = ExpMod)
  
  ############# Geometric Parameters ############
  Schematic_Eye$r_retina <- 11.06/23.01*Axial_Length; #retinal radius
  
  Schematic_Eye$a_retina <- 2*Schematic_Eye$r_retina; #Diameter of the retinal ellipsoid (major axis in the ellipse growth model)
  
  if (Schematic_Eye$ExpMod == 'Global'){
    Ell_Ratio <- 1; #This ratio makes the retinal ellipsoid a circle
  }
  
  if (Schematic_Eye$ExpMod == 'Elliptical'){
    pMod <- c(-0.0207, 1.4779); #Model for the ratio
    Ell_Ratio <- polyval(pMod, Schematic_Eye$a_retina);
  }
  
  Schematic_Eye$b_retina <- Ell_Ratio*Schematic_Eye$a_retina;
  Schematic_Eye$c_retina <- Axial_Length - Schematic_Eye$a_retina/2; #retinal centre
  
  Schematic_Eye$r_ant_L <- 10.0; #anterior lens radius
  Schematic_Eye$AL_Offset <- 3.60; #Anterior lens distance from vertex
  
  Schematic_Eye$r_post_L <- 6; #anterior lens radius
  Schematic_Eye$PL_Offset <- 7.375; #Posterior lens distance from vertex
  
  Schematic_Eye$corneal_ecc <- 0.5;#Eccentricity of corneal ellipse;
  Schematic_Eye$corneal_rad <- ((Schematic_Eye$C_Index - 1)*1000)/(K);#Radius of corneal ellipse;
  Schematic_Eye$alpha <- alpha
  
  Schematic_Eye$nodal_point <- nodal_point_calc(Schematic_Eye)
  
  return(Schematic_Eye)
}

##### Line-ellipse and circle intersection #######
lineEllipse <- function(a, b, h, k, alpha, p, q ){ 
  
  #Parameter of line
  m = -q/p;
  c = q;
  #Parameter for ellipse
  A = ( cos(alpha)^2 / a^2  + sin(alpha)^2 / b^2 );
  B = - 2 * cos (alpha) * sin(alpha) * (1/a^2 - 1/b^2);
  C = ( sin(alpha)^2 / a^2  + cos(alpha)^2 / b^2 );
  #Handle case with infinite slope
  if (q == Inf){
    x1 <- p;
    x2 <- p;
    M <- C;
    N <- B*p-(2*C*k+B*h);
    O <- A*p^2 - p*(2*A*h + k*B) + (A*h^2+B*h*k+C*k^2-1);
    determinant <- (N^2 - 4* M * O);
    y1  <- (-N + sqrt(determinant))/ (2*M);
    y2  <- (-N - sqrt(determinant))/ (2*M);
  }else{    
    M <- A + B*m + C* m^2;
    N <- B*c + 2*C*m*c - 2*A*h - k*B - m* (2*C*k+ B* h);
    O <- C*c^2 - c* (2*C*k + B*h) + A*h^2 + B*h*k + C*k^2 -1 ;
    
    determinant = (N^2 - 4* M * O);
    
    x1  <- (- N + sqrt(determinant))/ (2*M);
    x2  <- (- N - sqrt(determinant))/ (2*M);
    
    y1 = m*x1 + c;
    y2 = m*x2 + c;
  }
  
  cond1 <- ((((x1-h)*cos(alpha) + (y1-k)*sin(alpha))^2 /a^2 + ((x1-h)*sin(alpha) - (y1-k)*cos(alpha))^2/b^2) >= (1-1e-15)) & 
    ( (((x1-h)*cos(alpha) + (y1-k)*sin(alpha))^2/a^2 + ((x1-h)*sin(alpha) - (y1-k)*cos(alpha))^2/b^2) <= (1+1e-15));
  cond2 <-  ((((x2-h)*cos(alpha) + (y2-k)*sin(alpha))^2/a^2 + ((x2-h)*sin(alpha) - (y2-k)*cos(alpha))^2/b^2) >= (1-1e-15)) & 
    ( (((x2-h)*cos(alpha) + (y2-k)*sin(alpha))^2/a^2 + ((x2-h)*sin(alpha) - (y2-k)*cos(alpha))^2/b^2) <= (1+1e-15));
  
  # if (cond1 == 1 & cond2 == 0){
  #   x2 <- x1;
  #   y2 <- y1;
  # }else{
  #   if (cond1 == 0 & cond2 == 1){
  #     x1 <- x2;
  #     y1 <- y2;
  #   }else{ 
  #     if(cond1 == 0 & cond2 == 0){
  #       x1 <- NA;
  #       x2 <- NA;
  #       y1 <- NA;
  #       y2 <- NA;
  #     }
  #   }
  # }
  
  return(data.frame(x1,x2,y1,y2))
}

linecirc <- function(m, q, xc, yc, r){
  
  x1 <- -1 - xc
  x2 <- 1 - xc
  y1 <- (-m + q) - yc
  y2 <- (m + q) - yc
  
  dx <- x2 - x1
  dy <- y2 - y1
  dr <- sqrt(dx^2 + dy^2)
  D <- x1*y2 - x2*y1
  
  x1 <- (D*dy + sign(dy)*dx*sqrt(r^2*dr^2 - D^2))/dr^2
  x2 <- (D*dy - sign(dy)*dx*sqrt(r^2*dr^2 - D^2))/dr^2
  
  y1 <- (-D*dx + abs(dy)*sqrt(r^2*dr^2 - D^2))/dr^2
  y2 <- (-D*dx - abs(dy)*sqrt(r^2*dr^2 - D^2))/dr^2
  
  x <- unique(c(x1, x2)) + xc
  y <- unique(c(y1, y2)) + yc
  
  return(data.frame(x, y))
}

###### ray-tracing ######

###### ray-tracing ######
ray_tracing_eye <- function(slope, intercpt, Schematic_Eye){
  
  Ray_TraceX <- numeric(5);
  Ray_TraceY <- numeric(5);
  
  #Calculates the centers of circumferences that make up the anterior and
  #posterior face of the lens
  Cent_ALens <- Schematic_Eye$AL_Offset + Schematic_Eye$r_ant_L;
  Cent_PLens <- Schematic_Eye$PL_Offset - Schematic_Eye$r_post_L;
  
  #Distance between the centers
  d <- Cent_PLens - Cent_ALens;
  
  #Calculates where the two faces join (only for display, not important for
  #calculations)
  Lens_Edge_X <- -(d^2 - Schematic_Eye$r_ant_L^2 + Schematic_Eye$r_post_L^2)/(2*d);
  Lens_Edge_Y <- sqrt(Schematic_Eye$r_post_L^2 - Lens_Edge_X^2);
  Lens_Edge_X <- Lens_Edge_X + Cent_PLens;
  
  if (slope != 0 & ((abs(-intercpt/slope) > (Schematic_Eye$c_retina + Schematic_Eye$r_retina)))){
    #Does not calculate if the incident ray never crosses optical axis
    Exit_Angle <- NA;
    xoutV <- NA;
    youtV <- NA;
    Dist_on_Retina <- NA;
    Focal_Length <- NA;
    Lens_AngleIn <- NA;
  }else{
    if (slope == 0 & intercpt == 0){
      #If no angle between the ray and the optical axis, everything is 0
      Exit_Angle <- 0;
      xoutV <- 0;
      youtV <- 0;
      Dist_on_Retina <- 0;
      Focal_Length <- 0;
      Lens_AngleIn <- 0;
    }else{
      ##Calculations
      #Calculates the horizontal location of the corneal apex (ellipse)
      Cornea_offset <- -Schematic_Eye$corneal_rad/(1 - Schematic_Eye$corneal_ecc);
      
      
      
      #Calculates cornea
      y_cornea <- seq(-10,10,0.01);
      x_cornea <- -sqrt((Schematic_Eye$corneal_rad^2 - (1 - Schematic_Eye$corneal_ecc)*(y_cornea^2))/
                          (1 - Schematic_Eye$corneal_ecc)^2) - Cornea_offset;
      
      #Calcuates anterior an posterior faces of the lens
      CC <- useful::pol2cart(Schematic_Eye$r_ant_L, deg2rad(90:270));
      x_ant_L <- CC$x
      y_ant_L <- CC$y
      x_ant_L <- x_ant_L + Cent_ALens;
      
      CC <- useful::pol2cart(Schematic_Eye$r_post_L, deg2rad(-90:90));
      x_post_L <- CC$x
      y_post_L <- CC$y
      x_post_L = x_post_L + Cent_PLens;
      
      #Eliminates points beyond the point where faces join
      Keep_A <- sqrt((x_ant_L - (-Schematic_Eye$r_post_L + Schematic_Eye$PL_Offset))^2 + (y_ant_L)^2) <= (Schematic_Eye$r_post_L);
      Keep_P <- sqrt((x_post_L - (Schematic_Eye$AL_Offset + Schematic_Eye$r_ant_L))^2 + (y_post_L)^2) <= (Schematic_Eye$r_ant_L);
      
      x_ant_L <- x_ant_L[Keep_A];
      y_ant_L <- y_ant_L[Keep_A];
      x_post_L <- x_post_L[Keep_P];
      y_post_L <- y_post_L[Keep_P];
      
      #Calculates retinal sphere
      CC = useful::pol2cart(Schematic_Eye$r_retina, deg2rad(180:540));
      x_retina <- CC$x
      y_retina <- CC$y
      x_retina <- x_retina + Schematic_Eye$c_retina;
      Ell_Ratio <- Schematic_Eye$b_retina/Schematic_Eye$a_retina;
      y_retina <- Ell_Ratio*y_retina;
      
      #Cancels retina and corneal points beyond crossing
      #Keep = sqrt((x_cornea - Schematic_Eye.c_retina).^2 + (y_cornea).^2) >= (Schematic_Eye.r_retina);
      Keep <- !(x_cornea >  Schematic_Eye$AL_Offset);
      x_cornea <- x_cornea[Keep];
      y_cornea <- y_cornea[Keep];
      
      Keep <- (x_retina >  Schematic_Eye$AL_Offset);
      x_retina <- x_retina[Keep];
      y_retina <- y_retina[Keep];
      
      
      #################################
      
      #Calcuates a and b parameters of the corneal ellipse
      aEll <- -Cornea_offset;
      bEll <- sqrt(abs(Cornea_offset)*Schematic_Eye$corneal_rad);
      
      #Calculates the intersections between the incident ray and the
      #corneal ellipse (2 points)
      LEInt <- lineEllipse(aEll, bEll, -Cornea_offset, 0, 0, -intercpt/slope, intercpt);
      yC <- c(LEInt$y1, LEInt$y2);
      
      #Takes the most anterior intersection point
      yC <- yC[c(LEInt$x1, LEInt$x2) == min(c(LEInt$x1, LEInt$x2))];
      xC <- min(c(LEInt$x1, LEInt$x2));
      
      yC <- unique(yC);
      
      #Calculates the equation for the line tangent to the cornea at the
      #intersection
      xT <- (xC-1):(xC+1);
      qT <- ((aEll^2)*(bEll^2)/(aEll^2*yC));
      mT <-  -(bEll^2*(xC + Cornea_offset))/(aEll^2*yC);
      
      #Derives the perpendicular line
      mP <- -1/mT;
      qP <- yC - mP*(xC + Cornea_offset);
      
      #Uses the perpendicular line as the normal plane (atan(mP)) to 
      #apply Snell's law from air (n = 1) to cornea
      #Calculates a new line equation
      RefM_C <- Snell_law(atan(mP), slope, 1, Schematic_Eye$C_Index);
      RefQ_C <- yC - RefM_C*xC;
      
      #Calculates the intersections between the ray refracted from the 
      #cornea and the cicumference of the anterior face of the lens
      CC <- linecirc(RefM_C,RefQ_C,Cent_ALens,0,Schematic_Eye$r_ant_L);
      
      #Takes the most anterior intersection point
      youtAL <- CC$y[CC$x == min(CC$x)];
      xoutAL <- min(CC$x);
      
      #Kills the process if the refracted ray misses the lens entirely
      if (youtAL > Lens_Edge_Y){
        Exit_Angle <- NA;
        xoutV <- NA;
        youtV <- NA;
      }
      
      x_at_0 = -intercpt/slope;
      Ray_TraceX[1:2] <- c(-20, xC);
      Ray_TraceY[1:2] <- c(-20, xC)*slope + intercpt;
      
      
      ###################
      #Calculates the line equation of the perpendicular plane
      #through the crossing point at the anterior face of the lens. 
      #This is a circle now, so it is simply the radius through 
      #the incident point
      mPL1 <- (youtAL)/(xoutAL - Cent_ALens);
      qPL1 <- (youtAL) - mPL1*xoutAL;
      
      xT2 <- (xoutAL-1):(xoutAL+1);
      
      #Applies Snell's law between the aqueous and the
      #interior of the lens, through the anteior face of the lens.
      RefM_L1 <- Snell_law(atan(mPL1), RefM_C, Schematic_Eye$A_Index, Schematic_Eye$Lens_Index);
      RefQ_L1 <- youtAL - RefM_L1*xoutAL;
      
      #Calculates the intersections between the ray refracted from the
      #anterior face of the lens and the cicumference of the 
      #posterior face of the lens
      CC <- linecirc(RefM_L1,RefQ_L1,Cent_PLens,0,Schematic_Eye$r_post_L);
      youtPL <- CC$y[CC$x == max(CC$x)];
      xoutPL <- max(CC$x);
      
      Ray_TraceX[3:4] <- c(xoutAL, xoutPL);
      Ray_TraceY[3:4] <- c(xoutAL, xoutPL)*RefM_L1 + RefQ_L1;
      
      Lens_AngleIn <- rad2deg(atan(RefM_C));
      
      ###
      #Calculates the line equation of the perpendicular plane
      #through the crossing point at the posterior face of the lens.
      #This is a circle now, so it is simply the radius through
      #the incident point
      mPL2 <- (youtPL)/(xoutPL - Cent_PLens);
      qPL2 <- (youtPL) - mPL2*xoutPL;
      
      #Applies Snell's law between the interior of the lens and the
      #viterous, through the posterior face of the lens.
      RefM_L2 <- Snell_law(atan(mPL2), RefM_L1, Schematic_Eye$Lens_Index, Schematic_Eye$V_Index);
      RefQ_L2 <- youtPL - RefM_L2*xoutPL;
      
      #Calculates the intersections between the ray refracted from the
      #posterior face of the lens and the ellipsoid of the
      #retina
      #CC <- linecirc(RefM_L2,RefQ_L2,Schematic_Eye$c_retina,0,Schematic_Eye$r_retina);
      Retina_offset <- Schematic_Eye$c_retina;
      CC <- lineEllipse(Schematic_Eye$a_retina/2, Schematic_Eye$b_retina/2, Retina_offset, 0, 0, -RefQ_L2/RefM_L2, RefQ_L2);
      CC <- data.frame(x = c(CC$x1, CC$x2), y = c(CC$y1, CC$y2));
      
      youtV <- CC$y[CC$x == max(CC$x)];
      youtV <- youtV[1]
      xoutV <- max(CC$x);
      
      Exit_Angle <- rad2deg(atan(RefM_L2));
      
      #The distance on the retina is the arc on the retinal
      #circumference (angle = atan of the angular coefficient) - numerical for ellipsoid
      CircSect <- 2*asin(.5*sqrt(youtV^2 + (xoutV - (Schematic_Eye$c_retina + Schematic_Eye$r_retina))^2)/Schematic_Eye$r_retina);
      CircSect <- CircSect[1]
      
      if (Schematic_Eye$ExpMod == "Elliptical"){
        retinaArc <- useful::pol2cart(Schematic_Eye$r_retina, linspace(0,CircSect,100));
        Ell_Ratio <- Schematic_Eye$b_retina/Schematic_Eye$a_retina;
        x_retinaArc <- retinaArc$x
        y_retinaArc <- retinaArc$y
        y_retinaArc <- Ell_Ratio*y_retinaArc;
        Dist_on_Retina <- -sign(slope)*sum(sqrt((diff(y_retinaArc)^2 + diff(x_retinaArc)^2)));
      }
      
      if (Schematic_Eye$ExpMod == "Global"){
        Dist_on_Retina <- -sign(slope)*Schematic_Eye$r_retina*CircSect;
      }
      
      Focal_Length <- (-RefQ_L2/RefM_L2);
      
      Ray_TraceX[5] <- xoutV;
      Ray_TraceY[5] <- xoutV*RefM_L2 + RefQ_L2;
    }
  }
  
  Results <- list()
  Results[[1]] <- data.frame(Exit_Angle, Dist_on_Retina, xoutV, youtV, Focal_Length, Lens_AngleIn)
  Results[[2]] <- data.frame(Ray_TraceX, Ray_TraceY)
  
  return(Results)
}

###### calculate nodal point ######## NOT NEEDED
nodal_point_calc <- function(Schematic_Eye){
  
  error_nodal <- function(yoffset, Schematic_Eye, Incident_Angle){
    
    slope <- tan(deg2rad(Incident_Angle));
    
    # result <- tryCatch({
    #   Exit_Angle <- ray_tracing_eye(slope, yoffset, Schematic_Eye)
    # }, 
    # error = function(e){"error"})
    # 
    # if (result == "error"){
    #   Error <- 1000
    # }else{
    
    Exit_Angle <- ray_tracing_eye(slope, yoffset, Schematic_Eye);
    Exit_Angle <- Exit_Angle[[1]]$Exit_Angle
    
    if (is.na(Exit_Angle) | (yoffset <= 0)){
      Error <- 1000
    }else{
      Error <- abs(Incident_Angle - Exit_Angle)
    }
    #}
    
    return(Error);
  }
  
  
  TestAntSegment = (Schematic_Eye$C_Index == 1.336 & 
                      Schematic_Eye$A_Index == 1.336 & 
                      Schematic_Eye$V_Index == 1.336 & 
                      Schematic_Eye$Lens_Index == 1.43 & 
                      Schematic_Eye$AL_Offset == 3.6 &
                      Schematic_Eye$PL_Offset == 7.375 &
                      Schematic_Eye$r_ant_L == 10 &
                      Schematic_Eye$r_post_L == 6 &
                      Schematic_Eye$corneal_ecc == 0.5);
  
  if (TestAntSegment){
    return(6.9296)
  }else{
    Incident_Angle <- -seq(0.1,67,0.1);
    Nodal <- numeric(length(Incident_Angle));
    
    for (i in 1:length(Incident_Angle)){
      
      slope = tan(deg2rad(Incident_Angle[i]));
      StartOff = -slope*7;
      
      opt_Res <- nlminb(start = StartOff, objective = error_nodal, 
                        Incident_Angle = Incident_Angle[i],
                        Schematic_Eye = Schematic_Eye)
      
      Nodal[i] = -opt_Res$par/slope;
    }
    
    return(mean(Nodal))
  }
}

######## Visual degrees to mm calculation #####
Vis_deg2mm <- function(Schematic_Eye, LensRef = 0){
  #Estimate of the average nodal point (approximation, elliptical corneal
  #surface)
  if (!LensRef){
    NP <- Schematic_Eye$nodal_point #6.9296 #nodal_point_calc(Schematic_Eye);
  }else{
    NP <- Schematic_Eye.AL_Offset;
  }
  
  #Angles for numerical calculation
  Step <- 1;
  Angles <- seq(0,50,Step);
  
  #Reference point
  PointX <- NP;
  PointY <- 0;
  
  Dist_mm <- numeric(length(Angles));
  
  for (i in 1:length(Angles)){
    #Slope and intercept of the visual direction of the point in the visual 
    #field 
    Angle_Shift <- Angles[i];
    slope <- -tan(deg2rad(Angle_Shift))
    intercpt <- PointY - slope*PointX;
    
    #Ray tracing for each tested angle
    Ray_Trace <- ray_tracing_eye(slope, intercpt, Schematic_Eye);
    
    #Distances are calculated from the optic axis, always positive
    Dist_mm[i] <- Ray_Trace[[1]]$Dist_on_Retina
  }
  
  Area_ratio <- (diff(Dist_mm)^2)/(diff(Angles)^2)
  Angles <- Angles[1:(length(Angles)-1)] - Angles[1]
  Dist_mm = (Dist_mm[1:(length(Dist_mm)-1)]) - Dist_mm[1]
  
  return(data.frame(Angles, Dist_mm, Area_ratio))
}

############ 
FieldDeg2RetinaMM <- function(X, Y, Schematic_Eye){
  
  Conv <- Vis_deg2mm(Schematic_Eye)
  alpha <- Schematic_Eye$alpha
  
  #Assumes right eye, so shifts horizontally by alpha degrees nasally (negative)
  X <- X - alpha
  
  #Retinal arc corresponding to angle alpha
  Alpha_mm <- approx(Conv$Angles, Conv$Dist_mm, abs(alpha))$y
  
  Xmm <- sign(X)*approx(Conv$Angles, Conv$Dist_mm, abs(X))$y;
  Ymm <- sign(Y)*approx(Conv$Angles, Conv$Dist_mm, abs(Y))$y;
  
  Xmm <- Xmm + Alpha_mm
  
  CoordsMM <- data.frame(Xmm, Ymm)
  return(CoordsMM)
}

RetinaMM2FieldDeg <- function(X, Y, Schematic_Eye){
  
  Conv <- Vis_deg2mm(Schematic_Eye)
  alpha <- Schematic_Eye$alpha
  
  #Retinal arc corresponding to angle alpha
  Alpha_mm <- approx(Conv$Angles, Conv$Dist_mm, abs(alpha))$y
  
  #Assumes right eye, so shifts horizontally by alpha degrees nasally (negative)
  X <- X - Alpha_mm
  
  Xdeg <- sign(X)*approx(Conv$Dist_mm, Conv$Angles, abs(X))$y;
  Ydeg <- sign(Y)*approx(Conv$Dist_mm, Conv$Angles, abs(Y))$y;
  
  Xdeg <- Xdeg + alpha
  
  CoordsDeg <- data.frame(Xdeg, Ydeg)
  return(CoordsDeg)
}



########## THIS IS REQUIRED FOR DRASDO #################
load("Files/Disp_LUTNew.RData")

########## THIS IS REQUIRED FOR DRASDO #################

foveaDisc_Angle <- function(XFovea, YFovea, XONH, YONH){
  FovDisc_Angle <- atan((YFovea-YONH)/(XFovea-XONH))
  return(as.numeric(FovDisc_Angle))
}


drasdoFast <- function(X, Y, asdegree = 1, Ax_Length = (11.459/11.06*23.01), 
                       ONH = c(15, 0), ExpMod = 'Global',
                       SupPositive = TRUE, Eye = "Right", InvertD = FALSE){
  
  #Only one axial length calculated, left in for bacwards compatibility
  Ax_Length <- (11.459/11.06*23.01)
  
  #Drasdo displacement (from anatomical fovea). 
  #Default is assumed in retinal coordinates, right eye, in visual degrees, 
  #superior field indicated with postive values
  
  InvX <- 1
  InvY <- 1
  
  if (!SupPositive){
    Y <- -Y
    InvY <- -1
    ONH[2] <- -ONH[2]
  }
  if (tolower(Eye) == "left"){
    X <- -X
    InvX <- -1
    ONH[1] <- -ONH[1]
  }
  
  FovDisc_Angle <- foveaDisc_Angle(0, 0, ONH[1], ONH[2])
  
  Schematic_Eye <- schem_eye(Axial_Length = Ax_Length, ExpMod = ExpMod)
  
  #Rotate so that it matches Drasdo assumptions (ONH - Fovea angle is 0)
  R <- inv(rbind(c(cos(FovDisc_Angle), -sin(FovDisc_Angle)), c(sin(FovDisc_Angle), cos(FovDisc_Angle))))
  XY <- R%*%rbind(X, Y)
  X <- XY[1,]
  Y <- XY[2,]
  
  if (asdegree){
    CoordsMM <- FieldDeg2RetinaMM(X, Y, Schematic_Eye)
    Xmm <- CoordsMM$Xmm
    Ymm <- CoordsMM$Ymm
  }else{
    Xmm <- X;
    Ymm <- Y;
    
    CoordsDeg <- RetinaMM2FieldDeg(Xmm, Ymm, Schematic_Eye)
    X <- CoordsDeg$Xdeg
    Y <- CoordsDeg$Ydeg
  }
  
  if (InvertD){
    xdispMM <- bilinear(GX, GY, xlutInv, Ymm, Xmm)$z
    ydispMM <- bilinear(GX, GY, ylutInv, Ymm, Xmm)$z
  }else{
    xdispMM <- bilinear(GX, GY, xlut, Ymm, Xmm)$z
    ydispMM <- bilinear(GX, GY, ylut, Ymm, Xmm)$z
  }
  
  #Unchanged outside diplacement zone
  isOutside <- which(is.na(xdispMM) | (xdispMM == 0 & Xmm != 0) | (ydispMM == 0 & Ymm != 0))
  xdispMM[isOutside] <- Xmm[isOutside]
  ydispMM[isOutside] <- Ymm[isOutside]
  
  CoordsDeg <- RetinaMM2FieldDeg(xdispMM, ydispMM, Schematic_Eye)
  xdispDeg <- CoordsDeg$Xdeg
  ydispDeg <- CoordsDeg$Ydeg
  
  DistMM <- sqrt((xdispMM - Xmm)^2 + (ydispMM - Ymm)^2);
  DistDeg <- sqrt((xdispDeg - X)^2 + (ydispDeg - Y)^2);
  
  #Rerotate to match original input
  XY <- inv(R)%*%rbind(X, Y)
  X <- XY[1,]
  Y <- XY[2,]
  
  XY <- inv(R)%*%rbind(xdispDeg, ydispDeg)
  xdispDeg <- XY[1,]
  ydispDeg <- XY[2,]
  
  XY <- inv(R)%*%rbind(Xmm, Ymm)
  Xmm <- XY[1,]
  Ymm <- XY[2,]
  
  XY <- inv(R)%*%rbind(xdispMM, ydispMM)
  xdispMM <- XY[1,]
  ydispMM <- XY[2,]
  
  Results <- list()
  Results[[1]] <- data.frame(X = InvX*X, Y = InvY*Y, xdispDeg = InvX*xdispDeg, 
                             ydispDeg = InvY*ydispDeg, DistDeg)
  Results[[2]] <- data.frame(Xmm = InvX*Xmm, Ymm = InvY*Ymm, xdispMM = InvX*xdispMM, 
                             ydispMM = InvY*ydispMM, DistMM)
  
  return(Results)
}

G <- log10(c(pi*(0.1/2)^2, pi*(0.21/2)^2, pi*(0.43/2)^2, pi*(0.86/2)^2, pi*(1.72/2)^2));
patchDiam <- sqrt((10^(G[5]))/pi)*2

create_StimMaskEdge <- function(X, Y, patchDiam = 0.43){
  
  if((length(patchDiam) > 1) & (length(patchDiam) != length(X))){
    stop('Stimulus diameter must be either of length 1 or the same length as the grid coordinates')
  }else{
    
    D <- data.frame(X = X, Y = Y, patchDiam = patchDiam)
    
    xF <- numeric();
    yF <- numeric();
    Ind <- numeric();
    
    Seq <- seq(0,360,10);
    
    for (Points in 1:length(D$X)){
      Ind <- c(Ind, rep(Points,length(Seq)));
      AddX <- D$X[Points] + cos(deg2rad(Seq))*D$patchDiam[Points]/2;
      AddY <- D$Y[Points] + sin(deg2rad(Seq))*D$patchDiam[Points]/2;
      
      xF <- c(xF, AddX);
      yF <- c(yF, AddY);
    }
    
    return(data.frame(Ind, xF, yF))
  }
}

DistMaskAsEdges <- function(DistortedPatches, MaskEdges){
  DistortedDeg <- cbind(DistortedPatches[[1]][,3:4], MaskEdges$Ind)
  DistortedMM <- cbind(DistortedPatches[[2]][,3:4], MaskEdges$Ind)
  
  colnames(DistortedDeg) <- c("xF", "yF", "Ind")
  colnames(DistortedMM) <- c("xF", "yF", "Ind")
  
  List <- list()
  List[[1]] <- DistortedDeg
  List[[2]] <- DistortedMM
  
  return(List)
}

measure_StimMasks <- function(DistMasks, StructMap, tol = .5){
  Avrg <- numeric(length(unique(DistMasks$Ind)))*NA
  SumC <- numeric(length(unique(DistMasks$Ind)))*NA
  DArea <- numeric(length(unique(DistMasks$Ind)))*NA
  
  for (i in 1:length(unique(DistMasks$Ind))){
    Sel <- DistMasks$Ind == unique(DistMasks$Ind)[i]
    if (!is.na(sum(DistMasks$xF[Sel]))){
      w <- owin(poly = cbind(DistMasks$xF[Sel], DistMasks$yF[Sel]))
      SelectInside <- inside.owin(StructMap$X, StructMap$Y, w)
      NumTot <- sum(SelectInside)
      if (NumTot > 0){
        NumNas <- sum(is.na(StructMap$Z[SelectInside]))
        if ((NumNas/NumTot) <= 0.5){
          Avrg[i] <- mean(StructMap$Z[SelectInside], na.rm = TRUE)
          SumC[i] <- Avrg[i]*NumTot
          DArea[i] <- area.owin(w)
        }else{
          Avrg[i] <- NA
          SumC[i] <- NA
          DArea[i] <- NA
        }
        
      }
    }
  }
  
  MasksLog <- data.frame(Ind = unique(DistMasks$Ind), Avrg, SumC, DArea)
  
  return(MasksLog)
}


#Do not use
create_StimMasks <- function(MaskEdges, FoveaXY = c(ImSize/2, ImSize/2), ImSize = 768, FoVX = 30, FoVY = 30){
  
  DegToRad_X <- FoVX/ImSize
  DegToRad_Y <- FoVY/ImSize
  
  MaskXY <- list(x = (1:ImSize - FoveaXY[1])*DegToRad_X,
                 y = (1:ImSize - FoveaXY[2])*DegToRad_Y)
  
  MasksLog <- array(FALSE, c(ImSize, ImSize, length(unique(MaskEdges$Ind))))
  
  for (i in 1:length(unique(MaskEdges$Ind))){
    Sel <- MaskEdges$Ind == unique(MaskEdges$Ind)[i]
    w <- as.mask(owin(poly = cbind(MaskEdges$xF[Sel], MaskEdges$yF[Sel])), #eps = FoVX/ImSize, 
                 xy = MaskXY)
    MasksLog[,,i] <- w$m
    
  }
  
  return(MasksLog)
  
}
