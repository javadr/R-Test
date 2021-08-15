# The following code will demonstrate how the internal linear algebra functions of R are more efficient than user defined!
# It implements Gauss-Jordan algorithm in solving system of linear equations

# a very simple library for getting running's time of part of a code
#install.packages("tictoc")
library(tictoc)


# solving a system of linear equatinos with an internal function of R
solveLinearEquation <- function (A, b) {
  return (solve(A) %*% b)
}


# solving a system of linear equations via GaussJordan algorithm.
solveLinearEquationGaussJordan <- function (A, b) {
  if (det(A) == 0) {
    stop ("system is singular!")
  }
  
  n = dim(A)[1] # dimenstion of system; the number of linear equations
  linearSystem = cbind(A, b) # make new matrix joining A (equations) and b (answers)
  
  # iterates on equations top-down
  for (r in 1:(n - 1))
    # iterates on next equations
    for (br in (r + 1):n) {
      # in case getting zero, one of the following nonzero-rows must be exchanged with zero one!
      if (linearSystem[r, r] == 0) {
        nextNonZeroRow = r + which(linearSystem[, r][r + 1:n] != 0)
        linearSystem[c(r, nextNonZeroRow),] = linearSystem[c(nextNonZeroRow, r),]
      }
      if (linearSystem[br, r] != 0) {
        # alpha is a coeficient by that we can get to zero by multiplication and adding that row
        alpha = linearSystem[br, r] / linearSystem[r, r]
        linearSystem[br,] = linearSystem[br,] - linearSystem[r,] * alpha
      }
    }
  
  # iterates on rows bottom-up
  for (r in n:2)
    # iterates on previous rows
    for (ar in (r - 1):1)
      if (linearSystem[ar, r] != 0) {
        alpha = linearSystem[ar, r] / linearSystem[r, r]
        linearSystem[ar,] = linearSystem[ar,] - linearSystem[r,] * alpha
        
      }
  
  # making the final results by dividing the answers by the only number in equation
  for (i in 1:n)
    if (linearSystem[i, i] != 1) {
      linearSystem[i, ] = linearSystem[i, ] * (1 / linearSystem[i, i])
    }
  
  return (matrix(linearSystem[, (n + 1)]))
}

# A <- matrix(c(2, 4, 6, 4, 5, 6, 3, 1, -2), ncol = 3, byrow = TRUE)
# b <- c(18, 24, 4)

# A <- matrix(c(1, -1, -1, 1, 2, 0, 2, 0, 0, -1, -2, 0, 3, -3, -2, 4), ncol = 4, byrow = TRUE)
# b <- c(0, 8, -8, 7)

tic.clearlog()

for (i in 1:100) {
  n = sample(4:20, 1)
  A <- floor(matrix(runif(n * n) * 100, ncol = n))
  b <- floor(runif(n) * 100)
  
  # cat (sprintf("\niteration number %d\n", i))
  
  tic("R default solve function: ")
  resultR = solveLinearEquation(A, b)
  toc(log = TRUE, quiet = TRUE)
  
  tic("User defined function: ")
  resultUser =  solveLinearEquationGaussJordan(A, b)
  toc(log = TRUE, quiet = TRUE)
  
  if (!all(resultR - resultUser < .00001)) {
    cat ("THERE IS A MISMATCH!\n")
    print (resultR)
    print (resultUser)
  }
}

log.lst <- tic.log(format = FALSE)
tic.clearlog()

timings <- unlist(lapply(log.lst, function(x)
  x$toc - x$tic))

cat ("Averge Time Elapsed by R function is: ")
meanTimebyR = mean(timings[seq(1, by = 2, len = 100)])
print (meanTimebyR)

cat ("Averge Time Elapsed by User Defined function is: ")
meanTimebyUser = mean(timings[seq(2, by = 2, len = 100)])
print (meanTimebyUser)

cat ("\nFINAL RESULT:\n")
cat (
  sprintf(
    "User defined function is %.2f times slower than R function!",
    meanTimebyUser / meanTimebyR
  )
)
