#Required packages
Packages <- c('Matrix', 'spam', 'cPCG', 'matlib')
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
  LT = Matrix::chol(A, pivot = ..., cache = TRUE)
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

#### 2D #### 
t = system.time(function_cholsolv2d(A_2D(n)[[1]],A_2D(n)[[2]]), gcFirst = TRUE)

# Create sequences
n_2D <- (2^(2:7)+1)^2
nPlot_2D <- seq(1,
                16641,
                1)

# Save times in array
Times_2D <- c(0.08, 0.17, 0.44, 1.53, 5.64, 42.72)
# Create dataframe containing the saved times, and sequences of n_2D and its powers.
DF_2D <- data.frame(Time = Times_2D,
                    n = n_2D,
                    n2 = n_2D^2,
                    n3 = n_2D^3)
# Create second and third order polynomial models for the CPU times
Model_SO2D <- lm(Time ~ n + n2, data = DF_2D)
Model_TO2D <- lm(Time ~ n + n2 + n3, data = DF_2D)

# Calculate the fitted values from the models before
Fits_SO2D <- predict(Model_SO2D, list(n = nPlot_2D, n2 = nPlot_2D^2))
Fits_TO2D <- predict(Model_TO2D, list(n = nPlot_2D, n2 = nPlot_2D^2, n3 = nPlot_2D^3))

# Plot the second order polynomial fit and the CPU times
plot(n_2D, Times_2D, 
     main = 'Second order polynomial fit (2D)',
     ylab = 'CPU Time [s]',
     xlab = 'N')
lines(nPlot_2D, Fits_SO2D)

# Plot the third order polynomial fit and the CPU times
plot(n_2D, Times_2D, 
     main = 'Third order polynomial fit (2D)',
     ylab = 'CPU Time [s]',
     xlab = 'N')
lines(nPlot_2D, Fits_TO2D)

#### 3D ####
# Same concept as in 2D case
t = system.time(function_cholsolv3d(A_3D(n)[[1]],A_3D(n)[[2]]), gcFirst = TRUE)

n_3D <- (2^(2:5)+1)^3
nPlot_3D <- seq(1,
                35937,
                1)
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

#Incomplete Cholesky as basic iterative method
function_ICBIM = function(A, f, ...) {
  #initialization
  u = Matrix(0, nrow = nrow(f) , 1)
  I = diag(1, nrow = nrow(f), ncol = nrow(f))
  epsilon = 10^-10
  
  #calculate first resudual r
  r = f - A%*%u
  
  #calculate Incomplete Cholesky decomposition
  k = icc(as.matrix(A))
  M = k %*% t(k)
  Minv = inv(M)
  
  #calculate 2-norm of f
  normf = norm(f, type = '2')
  
  #calculate the convergence criterion
  convcrit = norm(r, type = '2')/normf
  iter = 0
  
  # convcrit_array <- c(convcrit)
  # iter_array <- c(iter)
  # r_array = c(r)
  
  #update solution and residual until convergence criterion is met
  while(epsilon <= convcrit){
    u = u + Minv %*% r
    r = (I - A %*% Minv) %*% r
    # r_array = c(r_array, r)
    convcrit = norm(r, type = '2')/normf
    # convcrit_array <- c(convcrit_array, convcrit)
    iter = iter + 1
    # iter_array <- c(iter_array, iter)
  }
  return(list(u = u, iter = iter_array, convcrit = convcrit_array, r = r_array))
}

# Plot for confirmation
ubim <- function_ICBIM(A_2D(9)[[1]], A_2D(9)[[2]])$u

n <- 9
x <- rep(seq(0, 1, by = 1/n),
         (n+1))
y <- x[order(x)]
uex <- function_uexact2d(x, y)

plot(ubim, main = 'IC as BIM (2D)', cex = 0, ylab = 'u', xlab = '')
legend(0, 1, legend=c("u IC as BIM", "u exact"),
       col=c("red", "blue"), lty=1:2, cex=0.8)
lines(ubim, col = 'red', lty = 1)
lines(uex, col = 'blue', lty = 2)


#### -------------------------- > 7 Log plot of condition against index -------------------------- ####

# 2D
Plot_List <- list()
j <- 0

for(i in c(5, 10, 15)){
  j <- j + 1
  ubim <- function_ICBIM(A_2D(i)[[1]], A_2D(i)[[2]])
  crit <- ubim$convcrit
  iter <- ubim$iter
  Plot_List[[j]] <- cbind(iter, crit)
}
plot(Plot_List[[3]], main = 'Convergence criterion as function of iteration', cex = 0, ylab = 'Ratio', xlab = 'm', log = 'y', ylim = c(10^(-11),1))
colours <- c('red', 'blue', 'black')
for(i in 1:length(Plot_List)){
  lines(Plot_List[[i]], col = colours[i], lty = i, lwd = 2)
}
abline(h = 10^(-10), col = 'green')
legend(100, 1, legend=c("N = 19", "N = 61", 'N = 127', 'Conv. Crit.'),
       col = c(colours, 'green'), lty = 1:3, cex=0.8)

# 3D
Plot_List <- list()
j <- 0

for(i in c(2, 4, 6)){
  j <- j + 1
  ubim <- function_ICBIM(A_3D(i)[[1]], A_3D(i)[[2]])
  crit <- ubim$convcrit
  iter <- ubim$iter
  Plot_List[[j]] <- cbind(iter, crit)
}
plot(Plot_List[[3]], main = 'Convergence criterion as function of iteration', cex = 0, ylab = 'Ratio', xlab = 'm', log = 'y', ylim = c(10^(-17),1))
colours <- c('red', 'blue', 'black')
for(i in 1:length(Plot_List)){
  lines(Plot_List[[i]], col = colours[i], lty = i, lwd = 2)
}
abline(h = 10^(-10), col = 'green')
legend(23, 1, legend=c("N = 3", "N = 15", 'N = 29', 'Conv. Crit.'),
       col = c(colours, 'green'), lty = 1:3, cex=0.8)

#### -------------------------- > 8 Residual reduction factor for last 5 iterations -------------------------- ####
# 2D
r_List <- list()
j <- 0

for(i in c(5, 10, 15)){
  j <- j + 1
  ubim <- function_ICBIM(A_2D(i)[[1]], A_2D(i)[[2]])
  r <- ubim$r
  r_List[[j]] <- r
}

rred_list <- list()
factor_list <- list()
for(i in 1:length(r_List)){
  rred_list[[i]] <- tail(sapply(r_List[[i]], norm, type = '2'), 5)
  factor_list[[i]] <- rred_list[[i]][2:5] / rred_list[[i]][1:4]
}

# 3D
r_List <- list()
j <- 0

for(i in c(2, 4, 6)){
  j <- j + 1
  ubim <- function_ICBIM(A_3D(i)[[1]], A_3D(i)[[2]])
  r <- ubim$r
  r_List[[j]] <- r
}

rred_list <- list()
factor_list <- list()
for(i in 1:length(r_List)){
  rred_list[[i]] <- tail(sapply(r_List[[i]], norm, type = '2'), 5)
  factor_list[[i]] <- rred_list[[i]][2:5] / rred_list[[i]][1:4]
}

#### -------------------------- > 9 CPU time for IC as BIM -------------------------- ####

#### 2D #### 

# Create sequences
n_2D <- ((1:16)+1)^2
nPlot_2D <- seq(1,
                289,
                1)

# Save times in array
Times_2D <- c(0.06, 0.08, 0.11, 0.14, 0.17, 0.25, 0.35, 0.42, 0.59, 0.82, 1.64, 2.83, 6.67, 10.47, 17.28, 25.11)
# Create dataframe containing the saved times, and sequences of n_2D and its powers.
DF_2D <- data.frame(Time = Times_2D,
                    n = n_2D,
                    n2 = n_2D^2,
                    n3 = n_2D^3)
# Create second and third order polynomial models for the CPU times
Model_SO2D <- lm(Time ~ n + n2, data = DF_2D)
Model_TO2D <- lm(Time ~ n + n2 + n3, data = DF_2D)

# Calculate the fitted values from the models before
Fits_SO2D <- predict(Model_SO2D, list(n = nPlot_2D, n2 = nPlot_2D^2))
Fits_TO2D <- predict(Model_TO2D, list(n = nPlot_2D, n2 = nPlot_2D^2, n3 = nPlot_2D^3))

# Plot the second order polynomial fit and the CPU times
plot(n_2D, Times_2D, 
     main = 'Second order polynomial fit (2D)',
     ylab = 'CPU Time [s]',
     xlab = 'N')
lines(nPlot_2D, Fits_SO2D)

# Plot the third order polynomial fit and the CPU times
plot(n_2D, Times_2D, 
     main = 'Third order polynomial fit (2D)',
     ylab = 'CPU Time [s]',
     xlab = 'N')
lines(nPlot_2D, Fits_TO2D)

#### 3D ####
# Same concept as in 2D case
t = system.time(function_cholsolv3d(A_3D(n)[[1]],A_3D(n)[[2]]), gcFirst = TRUE)

n_3D <- ((1:6)+1)^3
nPlot_3D <- seq(1,
                343,
                1)
Times_3D <- c(0.09, 0.14, 0.20, 0.52, 2.45, 5.84)
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


#### -------------------------- > 10 Ic as preconditioner -------------------------- ####

#Incomplete Cholesky as preconditioner on Conjugate Gradient method
function_ICCG = function(A, f, ...){
  #initialization
  u = Matrix(0, nrow = nrow(f) , 1)
  epsilon = 10^-10
  
  #calculate first resudual r
  r = f - A %*% u
  
  #calculate Incomplete Cholesky decomposition
  k = icc(as.matrix(A))
  M = k %*% t(k)
  Minv = inv(M)
  
  #calculate 2-norm of f
  normf = norm(f, type = '2')
  
  #calculate convergence criterion
  convcrit = norm(r, type = '2')/normf
  iter = 0
  
  convcrit_array <- c(convcrit)
  iter_array <- c(iter)
  
  #update values until convergence criterion is met
  while(epsilon <= convcrit){
    #preconditioning
    z = Minv %*% r
    
    if(iter == 0){
      p = z
    }
    else{
      beta = (t(r) %*% z)/(t(rold) %*% zold)
      p = z + beta[1] * p
    }
    
    #keep values r^(k-1) and z^(k-1) for the calculation of beta
    rold = r
    zold = z
    
    alpha = (t(r) %*% z) / (t(p) %*% A %*% p)
    u = u + alpha[1] * p
    r = r - alpha[1] * A %*% p
    
    iter = iter + 1
    convcrit = norm(r, type = '2')/normf
    
    convcrit_array <- c(convcrit_array, convcrit)
    iter_array <- c(iter_array, iter)
    
  } 
  
  return(list(u = u, iter = iter_array, convcrit = convcrit_array))
}

# Plot for confirmation
uiccg <- function_ICCG(A_3D(3)[[1]], A_3D(3)[[2]])$u

n <- 3
x <- rep(seq(0, 1, by = 1/n),
         (n+1)^2)
y <- rep(seq(0, 1, by = 1/n),
         n+1)
y <- rep(y[order(y)],
         n+1)
z <- rep(seq(0, 1, by = 1/n),
         (n+1)^2)
z <- z[order(z)]
uex <- function_uexact3d(x, y, z)

plot(uiccg, main = 'ICCG (3D)', cex = 0, ylab = 'u', xlab = '')
legend(0, 1, legend=c("u ICCG", "u exact"),
       col=c("red", "blue"), lty=1:2, cex=0.8)
lines(uiccg, col = 'red', lty = 1)
lines(uex, col = 'blue', lty = 2)

#### -------------------------- > 11 Log plot of condition against index, compare with 7 -------------------------- ####

# 2D
Plot_List <- list()
j <- 0

for(i in c(5, 10, 15)){
  j <- j + 1
  ubim <- function_ICCG(A_2D(i)[[1]], A_2D(i)[[2]])
  crit <- ubim$convcrit
  iter <- ubim$iter
  Plot_List[[j]] <- cbind(iter, crit)
}
plot(Plot_List[[3]], main = 'Convergence criterion as function of iteration', cex = 0, ylab = 'Ratio', xlab = 'm', log = 'y', ylim = c(10^(-12),1))
colours <- c('red', 'blue', 'black')
for(i in 1:length(Plot_List)){
  lines(Plot_List[[i]], col = colours[i], lty = i, lwd = 2)
}
abline(h = 10^(-10), col = 'green')
legend(16, 1, legend=c("N = 19", "N = 61", 'N = 127', 'Conv. Crit.'),
       col = c(colours, 'green'), lty = 1:3, cex=0.8)

# 3D
Plot_List <- list()
j <- 0

for(i in c(2, 4, 6)){
  j <- j + 1
  ubim <- function_ICCG(A_3D(i)[[1]], A_3D(i)[[2]])
  crit <- ubim$convcrit
  iter <- ubim$iter
  Plot_List[[j]] <- cbind(iter, crit)
}
plot(Plot_List[[3]], main = 'Convergence criterion as function of iteration', cex = 0, ylab = 'Ratio', xlab = 'm', log = 'y', ylim = c(10^(-17),1))
colours <- c('red', 'blue', 'black')
for(i in 1:length(Plot_List)){
  lines(Plot_List[[i]], col = colours[i], lty = i, lwd = 2)
}
abline(h = 10^(-10), col = 'green')
legend(9, 1, legend=c("N = 3", "N = 15", 'N = 29', 'Conv. Crit.'),
       col = c(colours, 'green'), lty = 1:3, cex=0.8)

#### -------------------------- > 12 CPU time for IC as preconditioner -------------------------- ####

#### 2D #### 

# Create sequences
n_2D <- ((1:16)+1)^2
nPlot_2D <- seq(1,
                289,
                1)

# Save times in array
Times_2D <- c(0.05, 0.07, 0.08, 0.10, 0.17, 0.17, 0.20, 0.27, 0.36, 0.52, 1.06, 3.18, 4.48, 7.86, 13.14, 24.03)
# Create dataframe containing the saved times, and sequences of n_2D and its powers.
DF_2D <- data.frame(Time = Times_2D,
                    n = n_2D,
                    n2 = n_2D^2,
                    n3 = n_2D^3)
# Create second and third order polynomial models for the CPU times
Model_SO2D <- lm(Time ~ n + n2, data = DF_2D)
Model_TO2D <- lm(Time ~ n + n2 + n3, data = DF_2D)

# Calculate the fitted values from the models before
Fits_SO2D <- predict(Model_SO2D, list(n = nPlot_2D, n2 = nPlot_2D^2))
Fits_TO2D <- predict(Model_TO2D, list(n = nPlot_2D, n2 = nPlot_2D^2, n3 = nPlot_2D^3))

# Plot the second order polynomial fit and the CPU times
plot(n_2D, Times_2D, 
     main = 'Second order polynomial fit (2D)',
     ylab = 'CPU Time [s]',
     xlab = 'N')
lines(nPlot_2D, Fits_SO2D)

# Plot the third order polynomial fit and the CPU times
plot(n_2D, Times_2D, 
     main = 'Third order polynomial fit (2D)',
     ylab = 'CPU Time [s]',
     xlab = 'N')
lines(nPlot_2D, Fits_TO2D)

#### 3D ####
# Same concept as in 2D case
t = system.time(function_cholsolv3d(A_3D(n)[[1]],A_3D(n)[[2]]), gcFirst = TRUE)

n_3D <- ((1:6)+1)^3
nPlot_3D <- seq(1,
                343,
                1)
Times_3D <- c(0.11, 0.14, 0.19, 0.42, 2.80, 6.28)
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







