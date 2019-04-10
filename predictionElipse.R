predictionEllipse <- function(mod, newdata, level = 0.95, ggplot = TRUE){
  # labels
  lev_lbl <- paste0(level * 100, "%")
  resps <- colnames(mod$coefficients)
  title <- paste(lev_lbl, "confidence ellipse for", resps[1], "and", resps[2])
  
  # prediction
  p <- predict(mod, newdata)
  
  # center of ellipse
  cent <- c(p[1,1],p[1,2])
  
  # shape of ellipse
  Z <- model.matrix(mod)
  Y <- mod$model[[1]]
  n <- nrow(Y)
  m <- ncol(Y)
  r <- ncol(Z) - 1
  S <- crossprod(resid(mod))/(n-r-1)
  
  # radius of circle generating the ellipse
  tt <- terms(mod)
  Terms <- delete.response(tt)
  mf <- model.frame(Terms, newdata, na.action = na.pass, 
                   xlev = mod$xlevels)
  z0 <- model.matrix(Terms, mf, contrasts.arg = mod$contrasts)
  rad <- sqrt((m*(n-r-1)/(n-r-m))*qf(level,m,n-r-m)*z0%*%solve(t(Z)%*%Z) %*% t(z0))
  
  # generate ellipse using ellipse function in car package
  ell_points <- car::ellipse(center = c(cent), shape = S, radius = c(rad), draw = FALSE)
  
  # ggplot2 plot
  if(ggplot){
    require(ggplot2, quietly = TRUE)
    ell_points_df <- as.data.frame(ell_points)
    ggplot(ell_points_df, aes(x, y)) +
      geom_path() +
      geom_point(aes(x = TOT, y = AMI), data = data.frame(p)) +
      labs(x = resps[1], y = resps[2], 
           title = title)
  } else {
    # base R plot
    plot(ell_points, type = "l", xlab = resps[1], ylab = resps[2], main = title)
    points(x = cent[1], y = cent[2])
  }
}
