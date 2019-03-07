f <- function(x){
  return( sin(x) - x )
}
fprime <- function(x){
  return( cos(x) - 1 )
}

secantMethod <- function(accuracy, p0, p1) {
  prev <- p0
  curr <- p1
  new <- 0
  pn <- c(prev, curr)
  while (abs(prev - curr) >= accuracy){
    new <- curr - (f(curr)*(curr-prev))/(f(curr)-f(prev))
    prev <- curr
    curr <- new
    pn <- c(pn, new)
    #print(curr)
  }
  return(pn)
}

newtonMethod <- function(accuracy, p0) {
  curr <- p0
  prev <-  curr*2 + accuracy*2  #make sure while loop runs at least once
  new <- 0
  pn <- curr
  while(abs(curr - prev) >= accuracy){
    new <- curr - f(curr)/fprime(curr)
    prev <- curr
    curr <- new
    pn <- c(pn, new)
    #print(curr)
  }
  return(pn)
}

S <- secantMethod(0.00001, pi/4, 3*pi/8)
N <- newtonMethod(0.00001, pi/4)
par(mfrow=c(1,2))   
x <- 0:(length(S)-1)
plot(x, S, type="p", col="red", xlab="n", ylab="pn", main="Secant Method pn")
x <- 1:length(N)
plot(x, N, type="p", col="green", xlab="n", ylab="pn", main="Newton Method pn")

