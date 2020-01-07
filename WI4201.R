#Required packages
Packages <- c('Matrix', 'spam', 'kernlab')
Packages_New <- Packages[!(Packages %in% installed.packages())]
install.packages(Packages_New)
invisible(lapply(Packages, library, character.only = TRUE))

#### -------------------------- > 1 Construct A^h and f -------------------------- ####

#### 2D ####
#Function f for 2D
function_FUN2 = function(x,y){
  return(-12*(x^2)*(y^5) - 20*(x^4)*(y^3))
}

#Boundary u for 2D
function_boundu2d = function(x,y){
  if(x == 0 || x == 1 || y == 0 || y == 1){
    return((x^4)*(y^5))
  }
  return()
}

#Construct matrix A for 1D problem
A_1D <- function(n){
  
  #initialization of h and values on band
  h <- 1/n
  singleband_value <- -1
  diagonal_value <- 2
  
  #Construct band matrix
  A <- 1/h^2 * bandSparse(n+1, n+1,
                          (-1):1,
                          list(rep(singleband_value, n), 
                               rep(diagonal_value, n+1), 
                               rep(singleband_value, n)))
  
  #Fix boundary values 
  A[1, 1] <- A[n+1, n+1] <- 1
  A[2, 1] <- A[1, 2] <- A[n,n+1] <- A[n+1, n] <- 0
  return(A)
}

#Construct matrix A and f for 2D problem
A_2D <- function(n){
  
  #initialize h 
  h <- 1/n
  
  #Construct identity with 0's at endpoints I_{n+1}
  I_Nminus1 <- Matrix(diag(c(0,
                             rep(1, n-1),
                             0),
                           nrow = n+1),
                      sparse = TRUE)
  
  #Apply tensor products with 1D solution, add them and update some values
  deel1 <- kronecker(A_1D(n), I_Nminus1)
  deel2 <- kronecker(I_Nminus1, A_1D(n))
  sum <- deel1 + deel2
  sum[1, 1] <- sum[n+1, n+1] <- sum[(n+1)^2-n, (n+1)^2-n] <- sum[(n+1)^2, (n+1)^2] <- 1
  
  # Initialise f
  f <- Matrix(matrix(0,
                     nrow = (n+1)^2,
                     ncol = 1),
              sparse = TRUE)
  
  #Fill in f for first and last row
  for(i in 1:(n+1)){
    f[i,1] = function_boundu2d((i-1)/n,0)
    f[n*(n+1) + i,1] = function_boundu2d((i-1)/n,1)
  }
  
  #Fill in f for interior rows, without eliminating boundary values
  for(i in 1:(n-1)){
    f[(i*(n+1) + 1),1] = function_boundu2d(0, (i/n))
    f[(i+1)*(n+1),1] = function_boundu2d(1, (i/n))
    for(j in 1:(n-1)){
      f[i*(n+1) + j+1, 1] = function_FUN2((j/n),(i/n))
    }
  }
  
  #Add boundary values to f
  for(i in 1:(n-1)){
    f[i*(n+1) + 2,1] = f[i*(n+1) + 2,1] + (1/h^2)*(f[i*(n+1) + 1,1]) 
    f[(i+1)*(n+1) -1,1] = f[(i+1)*(n+1) -1,1] + (1/h^2)*(f[(i+1)*(n+1),1]) 
    f[(n-1)*(n+1) +1 + i,1] = f[(n-1)*(n+1) +1 + i,1] + (1/h^2)*(f[n*(n+1) +1 + i,1]) 
    f[n+1 +1 + i,1] = f[n+1 +1 + i,1] + (1/h^2)*(f[1 + i,1]) 
  }
  
  return(list(sum,f))
}

#### 3D ####
#Function f for 3D
function_FUN4 = function(x,y,z){
  return(-12*(x^2)*(y^5)*(z^6) -20*(x^4)*(y^3)*(z^6) -30*(x^4)*(y^5)*(z^4))
}

#Boundary u for 3D
function_boundu3d = function(x,y,z){
  if(x == 0 || x == 1 || y == 0 || y == 1 || z == 0 || z == 1){
    return((x^4)*(y^5)*(z^6))
  }
  return()
}

A_3D <- function(n){
  
  #initialize h and indentity matrices
  h <- 1/n
  I_Nminus1 <- Matrix(diag(c(0,
                             rep(1, n-1),
                             0),
                           nrow = n+1),
                      sparse = TRUE)
  I2 <- Matrix(diag(1, nrow = n+1), sparse = TRUE)
  
  #Tensor product for -1 2 -1 band with -1 in x-direction
  deel1 <- kronecker(I_Nminus1, A_1D(n))
  deel2 <- kronecker(I_Nminus1, deel1)
  
  #Tensor product for -1 2 -1 band with -1 in y-direction
  deel3 <- kronecker(A_1D(n), I_Nminus1)
  deel4 <- kronecker(I_Nminus1, deel3)
  
  #Tensor product for -1 2 -1 band with -1 in z-direction
  deel5 <- kronecker(I_Nminus1, I_Nminus1)
  deel6 <- kronecker(A_1D(n), deel5)
  
  #Add tensor products
  sum <- deel2 + deel4 + deel6
  
  # Now set the value on the diagonals for boundary points to 1.
  # If the diagonal is 0, then set it to 0. ifelse() is vectorized.
  diag(sum) <- ifelse(diag(sum) == 0,
                      1,
                      diag(sum))
  
  #Initialize f
  f <- Matrix(matrix(0,
                     nrow = (n+1)^3,
                     ncol = 1),
              sparse = TRUE)
  
  #Fill in f for first and last plane
  for(i in 0:((n+1)^2 -1)){
    f[(i+1),1] = function_boundu3d(((i%%(n+1))/n),((i%/%(n+1))/n),0)
    f[((n+1)^3 - (n+1)^2 +1+i),1] = function_boundu3d(((i%%(n+1))/n),((i%/%(n+1))/n),1)
  }
  
  #Fill in f without eliminating boundary values
  for(i in 1:(n-1)){
    for(j in 1:(n+1)){
      f[i*(n+1)^2 +j,1] = function_boundu3d((j-1)/n,0,i/n)
      f[i*(n+1)^2 + n*(n+1) + j,1] = function_boundu3d((j-1)/n,1,i/n)
    }
    
    for(j in 1:(n-1)){
      f[(i*(n+1)^2 + (n+1)*j + 1),1] = function_boundu3d(0,(j/n),(i/n))
      f[(i*(n+1)^2 + (n+1)*(j+1)),1] = function_boundu3d(1,(j/n),(i/n))
      for(k in 1:(n-1)){
        f[(i*(n+1)^2 + (n+1)*j + k + 1),1] = function_FUN4((k/n),(j/n),(i/n))
      }
    }
  }
  
  #Add boundary values to f
  for(i in 2:n){
    for(j in 2:n){
      f[(i-1)*(n+1)^2 + n+1 + j,1] = f[(i-1)*(n+1)^2 + n+1 + j,1] + (1/h^2)*(f[(i-1)*(n+1)^2 + j,1])
      f[(i-1)*(n+1)^2 + (n-1)*(n+1) + j,1] = f[(i-1)*(n+1)^2 + (n-1)*(n+1) + j,1] + (1/h^2)*(f[(i-1)*(n+1)^2 + n*(n+1) + j,1]) 
      f[(i-1)*(n+1)^2 + (j-1)*(n+1) + 2,1] = f[(i-1)*(n+1)^2 + (j-1)*(n+1) + 2,1] + (1/h^2)*(f[(i-1)*(n+1)^2 + (j-1)*(n+1) + 1,1]) 
      f[(i-1)*(n+1)^2 + j*(n+1) - 1,1] = f[(i-1)*(n+1)^2 + j*(n+1) - 1,1] + (1/h^2)*(f[(i-1)*(n+1)^2 + j*(n+1),1]) 
      f[(n+1)^2 + (i-1)*(n+1) +j,1] = f[(n+1)^2 + (i-1)*(n+1) +j,1] + (1/h^2)*(f[(i-1)*(n+1) +j,1]) 
      f[(n-1)*(n+1)^2 + (i-1)*(n+1) +j,1] = f[(n-1)*(n+1)^2 + (i-1)*(n+1) +j,1] + (1/h^2)*(f[n*(n+1)^2 + (i-1)*(n+1) +j,1]) 
    }
  }
  
  return(list(sum,f))
}


#### -------------------------- > 2 Direct Solver -------------------------- ####
#### 2D ####

#Exact solution in 2D
function_uexact2d = function(x,y){
  return((x^4)*(y^5))
}

#Direct solver Cholesky decomposition 2D
function_cholsolv2d = function(A,f, ...){
  
  #Use chol function for Cholesky decomposition
  LT = Matrix::chol(A, pviot = ..., cache = TRUE)
  L = t(LT)
  
  #Forward solve Ly = f
  y = forwardsolve(L, f)
  
  #Backward solve Au = f
  u = backsolve(LT, y)
  
  return(u)
}

#Calculate accuracy of Cholesky decomposition in 2D
function_cholacc2d = function(u){
  
  #initialize vector as matrix for the exact solution
  uex = matrix(0, nrow = nrow(u), ncol = 1)
  n = sqrt(nrow(u)) -1
  
  #Fill in all exact solutions for u
  for(i in 0:(n)){
    uex[(1 + i*(n+1)),1] = function_uexact2d(0,(i/n))
    uex[(n+1 + i*(n+1)),1] = function_uexact2d(1,(i/n))
  }
  for(i in 2:n){
    uex[i,1] = function_uexact2d((i-1)/n,0)
    uex[n*(n+1) +i,1] = function_uexact2d((i-1)/n,1)
  }
  for(i in 1:(n-1)){
    for(j in 1:(n-1)){
      uex[(i*(n+1)+1+j),1] = function_uexact2d((j/n),(i/n))
    }
  }
  
  #Calculate the absolute norm
  accuracy = max(abs(u - uex))
  
  return(accuracy)
}

#### 3D ####
#Exact solution in 3D
function_uexact3d = function(x,y,z){
  return((x^4)*(y^5)*(z^6))
}

#Direct solver Cholesky decomposition 3D
function_cholsolv3d = function(A,f, ...){
  
  #initialize matrix A, f and the Cholesky decomposition
  A = as.matrix(as.data.frame(A))
  f = as.matrix(as.data.frame(f))
  LT = Matrix::chol(A, pivot = ..., cache = TRUE)
  L = t(LT)
  
  #Forward solve Ly = f
  y = forwardsolve(L, f)
  
  #Backward solve Au = f
  u = backsolve(LT, y)
  
  return(u)
}

#Calculate accuracy of Cholesky decomposition in 3D
function_cholacc3d = function(u){
  #Find n from the input vector
  n = round((nrow(u))^(1/3) -1)
  
  # Create sequence of numbers to denote the grid points and consider their ordering.
  x <- rep(seq(0, n, by = 1),
           (n+1)^2)
  y <- rep(seq(0, n, by = 1),
           n+1)
  y <- rep(y[order(y)],
           n+1)
  z <- rep(seq(0, n, by = 1),
           (n+1)^2)
  z <- z[order(z)]
  
  #Create a dataframe containing the sequences x y and z.
  Dataframe <- data.frame(cbind(z,y,x))
  
  #Apply function fuction_uexact3d to each row of the dataframe and store the result as a matrix.
  uex <- as.matrix(mapply(function_uexact3d, Dataframe$x/n, Dataframe$y/n, Dataframe$z/n))
  
  #Calculate the absolute norm
  accuracy = max(abs(u - uex))
  
  return(accuracy)
}

#### -------------------------- > 3 Plotting -------------------------- ####

# t = system.time(function_cholsolv2d(A_2D(n)[[1]],A_2D(n)[[2]]), gcFirst = TRUE)


#### 2D #### 
t = system.time(function_cholsolv2d(A_2D(n)[[1]],A_2D(n)[[2]]), gcFirst = TRUE)

n_2D <- 2^(2:7)
nPlot_2D <- seq(1,
                128,
                0.1)
Times_2D <- c(0.08, 0.17, 0.44, 1.53, 5.64, 42.72)
DF_2D <- data.frame(Time = Times_2D,
                    n = n_2D,
                    n2 = n_2D^2,
                    n3 = n_2D^3)
Model_SO2D <- lm(Time ~ n + n2, data = DF_2D)
Model_TO2D <- lm(Time ~ n + n2 + n3, data = DF_2D)

Fits_SO2D <- predict(Model_SO2D, list(n = nPlot_2D, n2 = nPlot_2D^2))
Fits_TO2D <- predict(Model_TO2D, list(n = nPlot_2D, n2 = nPlot_2D^2, n3 = nPlot_2D^3))

plot(n_2D, Times_2D, 
     main = 'Second order polynomial fit (2D)',
     ylab = 'CPU Time [s]',
     xlab = 'N')
lines(nPlot_2D, Fits_SO2D)

plot(n_2D, Times_2D, 
     main = 'Third order polynomial fit (2D)',
     ylab = 'CPU Time [s]',
     xlab = 'N')
lines(nPlot_2D, Fits_TO2D)

#### 3D ####
t = system.time(function_cholsolv3d(A_3D(n)[[1]],A_3D(n)[[2]]), gcFirst = TRUE)

n_3D <- 2^(2:5)
nPlot_3D <- seq(1,
                32,
                0.1)
Times_3D <- c(0.45, 1.58, 10.03, 290.27)
DF_3D <- data.frame(Time = Times_3D,
                    n = n_3D,
                    n2 = n_3D^2,
                    n3 = n_3D^3)

Model_SO3D <- lm(Time ~ n + n2, data = DF_3D)
Model_TO3D <- lm(Time ~ n + n2 + n3, data = DF_3D)

Fits_SO3D <- predict(Model_SO3D, list(n = nPlot_3D, n2 = nPlot_3D^2))
Fits_TO3D <- predict(Model_TO3D, list(n = nPlot_3D, n2 = nPlot_3D^2, n3 = nPlot_3D^3))

plot(n_3D, Times_3D, 
     main = 'Second order polynomial fit (3D)',
     ylab = 'CPU Time [s]',
     xlab = 'N')
lines(nPlot_3D, Fits_SO3D)

plot(n_3D, Times_3D, 
     main = 'Third order polynomial fit (3D)',
     ylab = 'CPU Time [s]',
     xlab = 'N')
lines(nPlot_3D, Fits_TO3D)


#### -------------------------- > 4 Fill-in ratio -------------------------- ####

n = 2^2
nzA = nnzero(A_2D(n)[[1]])
nzC = nnzero(t(chol(A_2D(n)[[1]])))
nzC/nzA


#### -------------------------- > 5 Reordering -------------------------- ####

n <- 2^7
(t = system.time(function_cholsolv2d(A_2D(n)[[1]], A_2D(n)[[2]], TRUE), gcFirst = TRUE))



#### -------------------------- > 6 IC as BIM -------------------------- ####

function_ICBIM2D = function(n) {
  u = Matrix(0, (n+1)^2 , 1)
  I = diag(1, (n+1)^2, 1)
  epsilon = 10^-10
  
  A = A_2D(n)[[1]]
  f = A_2D(n)[[2]]
  
  r = f - A%*%u
  
  ktrans = inchol(A)
  #M = transpose(ktrans) %*% ktrans
  
  #convcrit = norm(r)/norm(f)
  
  #while(convcrit <= epsilon){
  #u = u + inverse(M) %*% r
  #r = (I - A %*% inverse(M)) %*% r
  #convcrit = norm(r)/norm(f)
  #}
  
  return(ktrans)
}

function_ICBIM2D(3)

#### -------------------------- > 7 Log plot of condition against index -------------------------- ####

#### -------------------------- > 8 Residual reduction factor for last 5 iterations -------------------------- ####

#### -------------------------- > 9 CPU time for IC as BIM -------------------------- ####

#### -------------------------- > 10 Ic as preconditioner -------------------------- ####

function_ICCG2D = function(n){
  u = Matrix(0, (n+1)^2 , 1)
  A = A_2D(n)[[1]]
  f = A_2D(n)[[2]]
  epsilon = 10^-10
  iter = 1
  r = f - A%*%u
  Ktrans = ichol(A)
  M = transpose(ktrans) %*% ktrans
  convcrit = norm(r)/norm(f)
  
  while(convcrit <= epsilon){
    z = inverse(M)%*% r
    if(iter == 1){
      p = z
      iter = 2
    }
    else{
      beta = (transpose(r) %*% z)/(transpose(rold) %*% zold)
      p = z + beta * p
    }
    rold = r
    zold = z
    alpha = (transpose(r) %*% z) / (transpose(p) %*% A %*% p)
    u = u + inverse(M) %*% alpha * p
    r = r - alpha * A %*% p
    convcrit = norm(r)/norm(f)
  } 
  
  return(u)
}

#### -------------------------- > 11 Log plot of condition against index, compare with 7 -------------------------- ####

#### -------------------------- > 12 CPU time for IC as preconditioner -------------------------- ####












