## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  tidy.opts = list(width.cutoff = 70),
  tidy = TRUE
)


## ----eval=FALSE---------------------------------------------------------------
#  library(devtools)
#  install_github("thebrisklab/singR")
#  

## ---- eval=FALSE--------------------------------------------------------------
#  ld: warning: directory not found for option '-L/opt/R/arm64/gfortran/lib/gcc/aarch64-apple-darwin20.2.0/11.0.0'
#  ld: warning: directory not found for option '-L/opt/R/arm64/gfortran/lib'
#  ld: library not found for -lgfortran
#  clang: error: linker command failed with exit code 1 (use -v to see invocation)

## ---- eval=FALSE--------------------------------------------------------------
#  FLIBS =  -L/opt/R/arm64/gfortran/lib/gcc/aarch64-apple-darwin20.2.0/11.0.0 -L/opt/R/arm64/gfortran/lib -lgfortran -lemutls_w -lm

## ----eval=FALSE---------------------------------------------------------------
#  FLIBS =  -L/usr/local/gfortran/lib/gcc/aarch64-apple-darwin20.2.0/11.0.0 -L/usr/local/gfortran/lib -lgfortran -lm

## ----eval=FALSE, tidy=TRUE----------------------------------------------------
#  library(singR)
#  data(exampledata)
#  data <- exampledata
#  
#  lgrid = 33
#  par(mfrow = c(2,4))
#  # Components for X
#  image(matrix(data$sjX[1,], lgrid, lgrid), col = heat.colors(12),
#        xaxt = "n", yaxt = "n",main=expression("True S"["Jx"]*", 1"))
#  image(matrix(data$sjX[2,], lgrid, lgrid), col = heat.colors(12),
#        xaxt = "n", yaxt = "n",main=expression("True S"["Jx"]*", 2"))
#  image(matrix(data$siX[1,], lgrid, lgrid), col = heat.colors(12),
#        xaxt = "n", yaxt = "n",main=expression("True S"["Ix"]*", 1"))
#  image(matrix(data$siX[2,], lgrid, lgrid), col = heat.colors(12),
#        xaxt = "n", yaxt = "n",main=expression("True S"["Ix"]*", 2"))
#  
#  # Components for Y
#  image(vec2net(data$sjY[1,]), col = heat.colors(12), xaxt = "n", yaxt = "n",
#        main=expression("True S"["Jy"]*", 1"))
#  image(vec2net(data$sjY[2,]), col = heat.colors(12), xaxt = "n", yaxt = "n",
#        main=expression("True S"["Jy"]*", 2"))
#  image(vec2net(data$siY[1,]), col = heat.colors(12), xaxt = "n", yaxt = "n",
#        main=expression("True S"["Iy"]*", 1"))
#  image(vec2net(data$siY[2,]), col = heat.colors(12), xaxt = "n", yaxt = "n",
#        main=expression("True S"["Iy"]*", 2"))

## ----origin,echo=FALSE,out.width = "100%",fig.cap="True loadings in example 1."----
knitr::include_graphics("figs/Original.png",dpi = NA)

## ----eval=FALSE,tidy=TRUE-----------------------------------------------------
#  example1=singR(dX = data$dX,dY = data$dY,individual = T)
#  

## ----eval=FALSE, tidy=TRUE----------------------------------------------------
#  n.comp.X = NG_number(data$dX)
#  n.comp.Y = NG_number(data$dY)

## ----eval=FALSE, tidy=TRUE----------------------------------------------------
#  # JB on X
#  estX_JB = lngca(xData = data$dX, n.comp = n.comp.X, whiten = 'sqrtprec',
#                  restarts.pbyd = 20, distribution='JB')
#  Uxfull <- estX_JB$U
#  Mx_JB = est.M.ols(sData = estX_JB$S, xData = data$dX)
#  
#  # JB on Y
#  estY_JB = lngca(xData = data$dY, n.comp = n.comp.Y, whiten = 'sqrtprec',
#                  restarts.pbyd = 20, distribution='JB')
#  Uyfull <- estY_JB$U
#  My_JB = est.M.ols(sData = estY_JB$S, xData = data$dY)

## ----eval=FALSE, tidy=TRUE----------------------------------------------------
#  matchMxMy = greedymatch(scale(Mx_JB,scale = F), scale(My_JB,scale = F), Ux = Uxfull, Uy = Uyfull)
#  permJoint <- permTestJointRank(matchMxMy$Mx,matchMxMy$My)
#  joint_rank = permJoint$rj

## ----eval=FALSE, tidy=TRUE----------------------------------------------------
#  # Center X and Y
#  dX=data$dX
#  dY=data$dY
#  n = nrow(dX)
#  pX = ncol(dX)
#  pY = ncol(dY)
#  dXcentered <- dX - matrix(rowMeans(dX), n, pX, byrow = F)
#  dYcentered <- dY - matrix(rowMeans(dY), n, pY, byrow = F)
#  
#  # For X
#  # Scale rowwise
#  est.sigmaXA = tcrossprod(dXcentered)/(pX-1)
#  whitenerXA = est.sigmaXA%^%(-0.5)
#  xDataA = whitenerXA %*% dXcentered
#  invLx = est.sigmaXA%^%(0.5)
#  
#  # For Y
#  # Scale rowwise
#  est.sigmaYA = tcrossprod(dYcentered)/(pY-1)
#  whitenerYA = est.sigmaYA%^%(-0.5)
#  yDataA = whitenerYA %*% dYcentered
#  invLy = est.sigmaYA%^%(0.5)

## ----eval=FALSE, tidy=TRUE----------------------------------------------------
#  # Calculate the Sx and Sy.
#  Sx=matchMxMy$Ux[1:joint_rank,] %*% xDataA
#  Sy=matchMxMy$Uy[1:joint_rank,] %*% yDataA
#  
#  JBall = calculateJB(Sx)+calculateJB(Sy)
#  
#  # Penalty used in curvilinear algorithm:
#  rho = JBall/10

## ----eval=FALSE, tidy=TRUE----------------------------------------------------
#  # alpha=0.8 corresponds to JB weighting of skewness and kurtosis (can customize to use different weighting):
#  alpha = 0.8
#  #tolerance:
#  tol = 1e-10
#  
#  out <- curvilinear_c(invLx = invLx, invLy = invLy, xData = xDataA,
#                                   yData = yDataA, Ux = matchMxMy$Ux, Uy = matchMxMy$Uy,
#                                   rho = rho, tol = tol, alpha = alpha,
#                                   maxiter = 1500, rj = joint_rank)

## ----eval=FALSE, tidy=TRUE----------------------------------------------------
#  # Estimate Sx and Sy and true S matrix
#  Sjx = out$Ux[1:joint_rank, ] %*% xDataA
#  Six = out$Ux[(joint_rank+1):n.comp.X, ] %*% xDataA
#  Sjy = out$Uy[1:joint_rank, ] %*% yDataA
#  Siy = out$Uy[(joint_rank+1):n.comp.Y, ] %*% yDataA
#  
#  # Estimate Mj and true Mj
#  Mxjoint = tcrossprod(invLx, out$Ux[1:joint_rank, ])
#  Mxindiv = tcrossprod(invLx, out$Ux[(joint_rank+1):n.comp.X, ])
#  Myjoint = tcrossprod(invLy, out$Uy[1:joint_rank, ])
#  Myindiv = tcrossprod(invLy, out$Uy[(joint_rank+1):n.comp.Y, ])
#  
#  # signchange to keep all the S and M skewness positive
#  Sjx_sign = signchange(Sjx,Mxjoint)
#  Sjy_sign = signchange(Sjy,Myjoint)
#  Six_sign = signchange(Six,Mxindiv)
#  Siy_sign = signchange(Siy,Myindiv)
#  
#  Sjx = Sjx_sign$S
#  Sjy = Sjy_sign$S
#  Six = Six_sign$S
#  Siy = Siy_sign$S
#  
#  Mxjoint = Sjx_sign$M
#  Myjoint = Sjy_sign$M
#  Mxindiv = Six_sign$M
#  Myindiv = Siy_sign$M
#  
#  est.Mj=aveM(Mxjoint,Myjoint)
#  
#  trueMj <- data.frame(mj1=data$mj[,1],mj2=data$mj[,2],number=1:48)
#  SINGMj <- data.frame(mj1=est.Mj[,1],mj2=est.Mj[,2],number=1:48)
#  
#  

## ----eval=FALSE,echo=FALSE, tidy=TRUE-----------------------------------------
#  
#  lgrid = 33
#  par(mfrow = c(2,4))
#  
#  image(matrix(Sjx[1,], lgrid, lgrid), col = heat.colors(12), xaxt = "n", yaxt = "n",main=expression("Estimate S"["Jx"]*", 1"))
#  image(matrix(Sjx[2,], lgrid, lgrid), col = heat.colors(12), xaxt = "n", yaxt = "n",main=expression("Estimate S"["Jx"]*", 2"))
#  image(matrix(Six[1,], lgrid, lgrid), col = heat.colors(12), xaxt = "n", yaxt = "n",main=expression("Estimate S"["Ix"]*", 1"))
#  image(matrix(Six[2,], lgrid, lgrid), col = heat.colors(12), xaxt = "n", yaxt = "n",main=expression("Estimate S"["Ix"]*", 2"))
#  
#  image(vec2net(Sjy[1,]), col = heat.colors(12), xaxt = "n", yaxt = "n",main=expression("Estimate S"["Jy"]*", 1"))
#  image(vec2net(Sjy[2,]), col = heat.colors(12), xaxt = "n", yaxt = "n",main=expression("Estimate S"["Jy"]*", 2"))
#  image(vec2net(Siy[1,]), col = heat.colors(12), xaxt = "n", yaxt = "n",main=expression("Estimate S"["Iy"]*", 1"))
#  image(vec2net(Siy[2,]), col = heat.colors(12), xaxt = "n", yaxt = "n",main=expression("Estimate S"["Iy"]*", 2"))
#  
#  

## ----estiexample1,echo=FALSE,out.width = "100%",fig.cap="Estimated joint loadings in example 1."----
knitr::include_graphics("figs/Esti_example1.png",dpi = NA)

## ----eval=FALSE,echo=TRUE,tidy=TRUE-------------------------------------------
#  library(tidyverse)
#  library(ggpubr)
#  
#  t1 <- ggplot(data = trueMj)+
#    geom_point(mapping = aes(y=mj1,x=number))+
#    ggtitle(expression("True M"["J"]*", 1"))+
#    theme_bw()+
#    theme(panel.grid = element_blank())
#  
#  t2 <- ggplot(data = trueMj)+
#    geom_point(mapping = aes(y=mj2,x=number))+
#    ggtitle(expression("True M"["J"]*", 2"))+
#    theme_bw()+
#    theme(panel.grid = element_blank())
#  
#  #SING mj
#  
#  S1 <- ggplot(data = SINGMj)+
#    geom_point(mapping = aes(y=mj1,x=number))+
#    ggtitle(expression("Estimated M"["J"]*", 1"))+
#    theme_bw()+
#    theme(panel.grid = element_blank())
#  
#  S2 <- ggplot(data = SINGMj)+
#    geom_point(mapping = aes(y=mj2,x=number))+
#    ggtitle(expression("Estimated M"["J"]*", 2"))+
#    theme_bw()+
#    theme(panel.grid = element_blank())
#  
#  ggarrange(t1,t2,S1,S2,ncol = 2,nrow = 2)

## ----mjex1,echo=FALSE,out.width = "100%",fig.cap="Estimated joint subject scores in example 1."----
knitr::include_graphics("figs/MJ.png",dpi = NA)

