####################################################################
#
#       Chapter 2
#
####################################################################
############
# Figure 2.1
############
f <- function(x){sin(15*x)} # define our true function f(x)
xj <- seq(0.0,0.5,0.1) # define 6 equally spread points of known output x_j
D <- f(xj) # calculate vector D
# we then set Beta_0, sigma_u and theta accordingly
B0 = 0
sigu = 1.0
theta = 0.13
Ef <- B0 # E[f(x)]
varf <- sigu^2 # Var[f(x)]
ED <- rep(Ef,length(D)) # E[D]
covfD <- function(x){varf*exp((-abs(x-xj)^2)/theta^2)} # Cov[f(x),D]
# Calculate variance of D by the Cov[D,D]
varD <- c()
for (i in xj) {
  varD <- c(varD,covfD(i))
}
dim(varD) <- c(length(xj),length(xj)) # set it as a matrix
invvarD <- solve(varD) # Var[D]^-1
EDf <- function(x){Ef+covfD(x)%*%invvarD%*%(D-ED)} # Expectation of f(x) adjusted by D

# We now define the x range we wish to plot over, then calculate y=E_D[f(x)] over it
x <- seq(-0.03,0.53,0.0001)
y <- c()
for (i in x) {
  y <- c(y,c(EDf(i)))
}
# p represents the points of known output x_j's position on the plot
p <- c()
for (i in xj) {
  p <- c(p,c(EDf(i)))
}
# Define the limits and labels, then produce the plot
limx=c(0.0,0.51)
limy=c(-1.25,1.25)
labx="Rate Parameter Value x"
laby="Concentration of f(x)"
plot(x,y,type='l',col='blue',lwd=3,xlim=limx,ylim=limy, main = paste("Beta =", B0, ", Sigma =", sigu, ", Theta =", theta), xlab=labx,ylab=laby )

# Finally calculate the interval by defining Var_D[f(x)]
varDf <- function(x){varf-covfD(x)%*%invvarD%*%covfD(x)}
y2 <- c()
for (i in x) {
  y2 <- c(y2,c(varDf(i)))
}
# Calculate +/- 3*(st.dev.)
y3 <- y+3*(y2)^0.5
y4 <- y-3*(y2)^0.5
# Plot the lines
lines(x,y3,type='l',col='red',lwd=2)
lines(x,y4,type='l',col='red',lwd=2)
# And finally the points of known output
points(xj,p,lwd=6)

# Uncomment for second picture: Simply plot the line of f(x)
#lines(x,f(x),type='l',col='black',lwd=2)



############
# Figure 2.2
############
# We repeat the exact same process, but now with beta defined as a vector of the five beta_0 values
# And so we plot only the expectation line for these five values
# (previous figure must be run beforehand)
beta <- c(-11,-3.8,0,3.8,11)
Ef <- B0
ED <- rep(Ef,length(D))
EDf <- function(x){Ef+covfD(x)%*%invvarD%*%(D-ED)}
y <- c()
for (i in x) {
  y <- c(y,c(EDf(i)))
}
plot(x,y,type='l',col='blue',lwd=3,xlim=limx,ylim=limy, xlab=labx,ylab=laby )
#####
B0 = beta[1]
Ef <- B0
ED <- rep(Ef,length(D))
EDf <- function(x){Ef+covfD(x)%*%invvarD%*%(D-ED)}
y <- c()
for (i in x) {
  y <- c(y,c(EDf(i)))
}
lines(x,y,type='l',col='red',lwd=2)
#####
B0 = beta[2]
Ef <- B0
ED <- rep(Ef,length(D))
EDf <- function(x){Ef+covfD(x)%*%invvarD%*%(D-ED)}
y <- c()
for (i in x) {
  y <- c(y,c(EDf(i)))
}
lines(x,y,type='l',col='green',lwd=2)
#####
B0 = beta[4]
Ef <- B0
ED <- rep(Ef,length(D))
EDf <- function(x){Ef+covfD(x)%*%invvarD%*%(D-ED)}
y <- c()
for (i in x) {
  y <- c(y,c(EDf(i)))
}
lines(x,y,type='l',col='yellow',lwd=2)
#####
B0 = beta[5]
Ef <- B0
ED <- rep(Ef,length(D))
EDf <- function(x){Ef+covfD(x)%*%invvarD%*%(D-ED)}
y <- c()
for (i in x) {
  y <- c(y,c(EDf(i)))
}
lines(x,y,type='l',col='brown',lwd=2)
#####
p <- c()
for (i in xj) {
  p <- c(p,c(EDf(i)))
}
points(xj,p,lwd=6)

#lines(x,z,type='l',col='black',lwd=3)
legend(0.0, -0.4, legend=beta,
       col=c("red","green","blue","yellow","brown"), lty=1, cex=0.8)



############
# Figure 2.3
############
# Similarly to 2.2, we define sigs a vector of the different sigma_u values
# But now only plot the interval lines
# (Figure 2.1 must be run beforehand)
sigs <- c(0.2,0.6,1.0,1.2,1.4)
f <- function(x){sin(15*x)}
xj <- seq(0.0,0.5,0.1)
D <- f(xj)
B0 = 0
sigu = sigs[1]
theta = 0.13
Ef <- B0
varf <- sigu^2
ED <- rep(Ef,length(D))
covfD <- function(x){varf*exp((-abs(x-xj)^2)/theta^2)}
varD <- c()
for (i in xj) {
  varD <- c(varD,covfD(i))
}
dim(varD) <- c(length(xj),length(xj))
invvarD <- solve(varD)
EDf <- function(x){Ef+covfD(x)%*%invvarD%*%(D-ED)}
y <- c()
for (i in x) {
  y <- c(y,c(EDf(i)))
}
p <- c()
for (i in xj) {
  p <- c(p,c(EDf(i)))
}
plot(x,y,type='l',col='white',xlim=limx,ylim=limy,lwd=2,xlab=labx,ylab=laby, main="Changing Values of Sigma")

varDf <- function(x){varf-covfD(x)%*%invvarD%*%covfD(x)}
y2 <- c()
for (i in x) {
  y2 <- c(y2,c(varDf(i)))
}
y3 <- y+3*(y2)^0.5
y4 <- y-3*(y2)^0.5
lines(x,y3,type='l',col='red',lwd=2)
lines(x,y4,type='l',col='red',lwd=2)

#####
sigu = sigs[2]
varf <- sigu^2
covfD <- function(x){varf*exp((-abs(x-xj)^2)/theta^2)}
varD <- c()
for (i in xj) {
  varD <- c(varD,covfD(i))
}
dim(varD) <- c(length(xj),length(xj))
invvarD <- solve(varD)
varDf <- function(x){varf-covfD(x)%*%invvarD%*%covfD(x)}
y2 <- c()
for (i in x) {
  y2 <- c(y2,c(varDf(i)))
}
y3 <- y+3*(y2)^0.5
y4 <- y-3*(y2)^0.5
lines(x,y3,type='l',col='green',lwd=2)
lines(x,y4,type='l',col='green',lwd=2)
#####
sigu = sigs[3]
varf <- sigu^2
covfD <- function(x){varf*exp((-abs(x-xj)^2)/theta^2)}
varD <- c()
for (i in xj) {
  varD <- c(varD,covfD(i))
}
dim(varD) <- c(length(xj),length(xj))
invvarD <- solve(varD)
varDf <- function(x){varf-covfD(x)%*%invvarD%*%covfD(x)}
y2 <- c()
for (i in x) {
  y2 <- c(y2,c(varDf(i)))
}
y3 <- y+3*(y2)^0.5
y4 <- y-3*(y2)^0.5
lines(x,y3,type='l',col='blue',lwd=2)
lines(x,y4,type='l',col='blue',lwd=2)
#####
sigu = sigs[4]
varf <- sigu^2
covfD <- function(x){varf*exp((-abs(x-xj)^2)/theta^2)}
varD <- c()
for (i in xj) {
  varD <- c(varD,covfD(i))
}
dim(varD) <- c(length(xj),length(xj))
invvarD <- solve(varD)
varDf <- function(x){varf-covfD(x)%*%invvarD%*%covfD(x)}
y2 <- c()
for (i in x) {
  y2 <- c(y2,c(varDf(i)))
}
y3 <- y+3*(y2)^0.5
y4 <- y-3*(y2)^0.5
lines(x,y3,type='l',col='yellow',lwd=2)
lines(x,y4,type='l',col='yellow',lwd=2)
#####
sigu = sigs[5]
varf <- sigu^2
covfD <- function(x){varf*exp((-abs(x-xj)^2)/theta^2)}
varD <- c()
for (i in xj) {
  varD <- c(varD,covfD(i))
}
dim(varD) <- c(length(xj),length(xj))
invvarD <- solve(varD)
varDf <- function(x){varf-covfD(x)%*%invvarD%*%covfD(x)}
y2 <- c()
for (i in x) {
  y2 <- c(y2,c(varDf(i)))
}
y3 <- y+3*(y2)^0.5
y4 <- y-3*(y2)^0.5
lines(x,y3,type='l',col='brown',lwd=2)
lines(x,y4,type='l',col='brown',lwd=2)
#####
points(xj,p,lwd=6)

legend(0.0,-0.5, legend=sigs,
       col=c("red","green","blue","yellow","brown"), lty=1, cex=0.8)



############
# Figure 2.4
############
# These plots are done exactly as 2.2/2.3, except with vector defined
thetas <- c(0.04,0.11,0.13,0.15,0.55)
# And the expectations and standard deviations calculated accordingly.



############
# Figure 2.5
############
# Produced with code of figure 2.1, with
B0 = 0 # then B0 = 1  
sigu = 1.0
theta = 0.01



############
# Figure 2.6
############
# First run figure 2.1
# set the z value we want to match then plot its line
z=-0.75 
zs<-rep(z,length(x))
lines(x,zs,lwd=2)
###
varee=(0.01*2)^2 # Take Var[epsilon]=0.0001
# Define and calculate implausibility
I2<- function(x){((EDf(x)-z)^2)/(varDf(x)+varee)}
y<-c()
for (i in x) {
  y <- c(y,c(I2(i))^0.5)
}
# Calculate which x values have implausibilty < 3
x2<-c()
for (i in x) { if (c(I2(i))^0.5 < 3){
  x2 <- c(x2,i)
}}
# And plot these along the bottom
y2<-rep(-1.34,length(x2))
points(cbind(x2,y2),col='green')
###
plot(x,y,xlim=limx,type='l',lwd=2,ylim=c(0,10),xlab=labx,ylab="Implausibility") # Plot Implausibility
ds<-rep(3,length(x)) # Plot Implausibility=3
lines(x,ds,lwd=2,col='red')
###
# Pick points to be added by finding the middle of each interval, then add them to a new vector xj2
interval1<-x2[which(0.2<x2 & x2<0.3)]
median(interval1) #(mean gets the same result!)
#0.26145
interval2<-x2[which(0.3<x2 & x2<0.4)]
median(interval2)
#0.3518
added<-c(0.26145,0.3518)
xj2<-c()
for(i in seq(1,length(xj),1)){
  if(xj[i]>added[1]){
    xj2<-c(xj[1:i-1],added[1],added[2],xj[i:length(xj)])
    stop()
  }
}
xj<-xj2 # Set xj as xj2
# Now the second two plots can be made, running the code of figure 2.1 WITHOUT redefining xj (starting at D<-f(xj))
# then replotting implausibility as above





####################################################################
#
#       Chapter 3
#
####################################################################
# For this chapter we create a function that allows us to plot all different plots at once or choose which ones
# to plot. Some terms however must be defined outside the function. 
# We show a worked example based on figure 3.1

# First define all parameter values
B0 = 0.25
sigu = 1.5
theta = 0.8
# Then we set values along the x and y axes, and form them into a grid
x <- seq(1.2,3.4,0.01)
y <- seq(1.2,3.4,0.01)
xy <- expand.grid(x,y)
# Define our toy function
f <- function(x){cos(0.5*x[,1])+sin(x[,1]*x[,2])} 
# Create a grid of points of known output
xjs <- seq(1.2,3.4,2.2/3)
yjs <- seq(1.2,3.4,2.2/3)
js <- as.matrix(expand.grid(xjs,yjs))

# Now we must create a function that uses these to create suitable plots:
# We call this fp2, and it takes TRUE/FALSE values for the first 5 inputs to decide whether to plot the true plot,
# the expectation plot, the standard deviation plot, the implausibility plot and the diagnostic plot.
# Levelss defines the contour levels of the true and expectation plots.
# ssave take TRUE/FALSE and defines whether to automatically save the plots or not
# nahme decides the name under which to save any plots
fp2<-function(actualplot,expectation,SD,implausibility,SK,Levelss,ssave,nahme){
  xy2<-cbind(xy[,1],xy[,2]) # Creates a matrix where each row represents a coordinate, defining every coordinate on the plane.
  if(actualplot==TRUE){
    z<-f(xy2) # creates z as all the values of f(x) at each coordinate
    dim(z) <- c(length(x),length(y)) # defines z as a matrix with values matching coordinate position
    
    if(ssave==TRUE){pdf(file=file.path("C:","Users","tcvin","OneDrive","Work","DIS","pics","core",paste(nahme,"a",".pdf",sep="")))} 
    # saves plot as pdf with a added for actualplot, obviously path must be changed for individual user
    filled.contour(x,y, z, color = rainbow, # plots true function as a rainbow contour
                   plot.title = {title(main = "Title",
                                       xlab = expression('x'[1]), ylab = expression('x'[2]))
                     points(js[,1],js[,2],lwd=6)
                   },
                   levels = Levelss,
                   key.title = title(main = "key"))
    if(ssave==TRUE){dev.off()}
    
  }
  D <- f(js) # Defines vector of points of known output, D
  Ef <- B0 # E[f(x)]
  varf <- sigu^2 # Var[f(x)]
  ED <- rep(Ef,length(D)) # E[D]
  # We create a Cov[f(x),D] function that takes input of matrix x with two rows containing multiple coordinates with which to use. It then uses the
  # dist function to calculate the distances between these coordinates and points js, carefully ignoring the values between the coordinates themselves
  # and points js themselves. It then transforms these according the the necessary equation. Each column represents the Cov[f(x),D] vector for each coordinate.
  covfD <- function(x){ 
    varf*exp(-(as.matrix(dist(rbind(x,js))))[(1+dim(x)[1]):(dim(x)[1]+dim(js)[1]),1:dim(x)[1]]^2/theta^2)
  }
  varD <- covfD(js) # Var[D]=Cov[D,D]
  invvarD <- solve(varD) # Var[D]^-1
  #EDf <- function(x){Ef+covfD(x)%*%invvarD%*%(D-ED)}
  # We create a E_D[f(x)] function that essentially calculates Ef+covfD(x)%*%invvarD%*%(D-ED) but more efficiently for large matrices
  EDf <- function(x){
    A<-covfD(x)
    C<-rep(D-ED,dim(x)[1])
    dim(C)<-c(length(D-ED),dim(x)[1])
    B<-invvarD%*%C
    Ef+colSums(A * B)
  }
  #Similarly for Var_D[f(x)] = varf-covfD(x)%*%invvarD%*%c(covfD(x))
  varDf <- function(x){
    A<-covfD(x)
    B<-invvarD%*%A
    varf-colSums(A * B)
  }
  if(expectation==TRUE){
    # Here we essentially calculate z<-EDf[xy2], but as discussed split the calculation into sets of 100 coordinates for efficiency.
    # As before we then produce a contour plot of this z.
    z<-c()
    for (i in 0:floor((dim(xy2)[1]/100)-1)) {
      z<-c(z,c(EDf(xy2[(i*100):(i*100+99),])))
    }
    q<-dim(xy2)[1]-100*floor((dim(xy2)[1]/100)-1)-100
    z<-c(z,c(EDf(xy2[(dim(xy2)[1]-q):(dim(xy2)[1]),])))
    dim(z) <- c(length(x),length(y))
    if(ssave==TRUE){pdf(file=file.path("C:","Users","tcvin","OneDrive","Work","DIS","pics","core",paste(nahme,"a",".pdf",sep="")))}
    filled.contour(x,y, z, color = rainbow,
                   plot.title = {title(main = "Title",
                                       xlab = expression('x'[1]), ylab = expression('x'[2]))
                     points(js[,1],js[,2],lwd=6)
                   },
                   levels = Levelss,
                   key.title = title(main = "key"))
    if(ssave==TRUE){dev.off()}
  }
  if(SD==TRUE){
    # Similarly for standard deviation
    z<-c()
    for (i in 0:floor((dim(xy2)[1]/100)-1)) {
      z<-c(z,c(varDf(xy2[(i*100):(i*100+99),])^0.5))
    }
    q<-dim(xy2)[1]-100*floor((dim(xy2)[1]/100)-1)-100
    z<-c(z,c(varDf(xy2[(dim(xy2)[1]-q):(dim(xy2)[1]),])^0.5))
    dim(z) <- c(length(x),length(y))    
    if(ssave==TRUE){pdf(file=file.path("C:","Users","tcvin","OneDrive","Work","DIS","pics","core",paste(nahme,"sd",".pdf",sep="")))}
    filled.contour(x,y, z, color = terrain.colors,
                   plot.title = {title(main = "Title",
                                       xlab = expression('x'[1]), ylab = expression('x'[2]))
                     points(js[,1],js[,2],lwd=6)
                   },
                   levels = seq(0,1,0.05),
                   key.title = title(main = "key"))
    if(ssave==TRUE){dev.off()}
  }
  if(implausibility==TRUE){
    # For implausibilty we must also define Var(epsilon)=varee and the value we are trying to match = valuu outside the function.
    # We then define a function for Implausibility^2 and plot these values over the plane as before.
    I2<- function(x){((EDf(x)-valuu)^2)/(varDf(x)+varee)}
    z<-c()
    for (i in 0:floor((dim(xy2)[1]/100)-1)) {
      z<-c(z,c(I2(xy2[(i*100):(i*100+99),])^0.5))
    }
    q<-dim(xy2)[1]-100*floor((dim(xy2)[1]/100)-1)-100
    z<-c(z,c(I2(xy2[(dim(xy2)[1]-q):(dim(xy2)[1]),])^0.5))
    dim(z) <- c(length(x),length(y))
    
    if(ssave==TRUE){pdf(file=file.path("C:","Users","tcvin","OneDrive","Work","DIS","pics","core",paste(nahme,"i",".pdf",sep="")))}
    filled.contour(x,y, z, color = colorRampPalette(c('cyan','brown','red'),1.6),
                   levels = c(0,1,2,3,7,12,18,25,33,42,52),
                   plot.title = {title(main = "Title",
                                       xlab = expression('x'[1]), ylab = expression('x'[2]))
                     points(js[,1],js[,2],lwd=6)
                   },
                   key.title = title(main = "key"))
    if(ssave==TRUE){dev.off()}
  }
  if(SK==TRUE){
    # We define SC as a function that calculates the diagnostic value for each coordinate, then again plot this over the plane
    SC<- function(x){(EDf(x)-f(x))/(varDf(x))^0.5}
    z<-c()
    for (i in 0:floor((dim(xy2)[1]/100)-1)) {
      z<-c(z,c(SC(xy2[(i*100):(i*100+99),])))
    }
    q<-dim(xy2)[1]-100*floor((dim(xy2)[1]/100)-1)-100
    z<-c(z,c(SC(xy2[(dim(xy2)[1]-q):(dim(xy2)[1]),])))
    dim(z) <- c(length(x),length(y))
    
    if(ssave==TRUE){pdf(file=file.path("C:","Users","tcvin","OneDrive","Work","DIS","pics","core",paste(nahme,"sk",".pdf",sep="")))}
    
    filled.contour(x,y, z, color = colorRampPalette(c('red','chartreuse','navy'),1),
                   levels = seq(-2,2,0.1),
                   plot.title = {title(main = "Title",
                                       xlab = expression('x'[1]), ylab = expression('x'[2]))
                     points(js[,1],js[,2],lwd=6)
                   },
                   key.title = title(main = "key"))
    if(ssave==TRUE){dev.off()}
  }
}

ss<-FALSE
# we set this so that for all plots we do not save


####################
# Figures 3.1 & 3.2
####################
x <- seq(1.2,3.4,0.01)
y <- seq(1.2,3.4,0.01)
xy <- expand.grid(x,y)
f <- function(x){cos(0.5*x[,1])+sin(x[,1]*x[,2])}
xjs <- seq(1.2,3.4,2.2/3)
yjs <- seq(1.2,3.4,2.2/3)
js <- as.matrix(expand.grid(xjs,yjs))
B0 = 0.25
sigu = 1.5
theta = 0.8
fp2(actual=TRUE,expectation=TRUE,SD=TRUE,implausibility=FALSE,SK=TRUE,Levelss=seq(-1.2,2,0.2),ssave=ss,nahme = 'cossing')



####################
# Figures 3.1 & 3.2
####################
x <- seq(0,0.5,length=81)
y <- seq(-1,1,0.025)
xy <- expand.grid(x,y)
f <- function(x){x[,2]*sin(16*x[,1])}
B0 = 0
sigu = 1.5
theta = 0.8
library(lhs)
rl <- maximinLHS(16,2) # Because this was done randomly the individual points and therefore plot will be different
js <- cbind(rl[,1]*0.5,rl[,2]*2-1) #shift/stretch the [0,1]*[0,1] rl to fit the plane
fp2(actual=TRUE,expectation=TRUE,SD=FALSE,implausibility=FALSE,SK=FALSE,Levelss=seq(-1.1,1.1,0.1),ssave=ss,nahme = 'badysin')

# transform:
x <- seq(0,2,0.025)
y <- seq(-1,1,0.025)
xy <- expand.grid(x,y)
f <- function(x){x[,2]*sin(4*x[,1])}
js <- cbind(rl[,1]*2,rl[,2]*2-1) #transform the same rl to fit the transformed plane
fp2(actual=TRUE,expectation=TRUE,SD=FALSE,implausibility=FALSE,SK=FALSE,Levelss=seq(-1.1,1.1,0.1),ssave=ss,nahme = 'goodysin')
# This still depends on whether you get a "good" maximinLHS - hence why we wrote our own function. 



##################
# Figure 3.4
##################
# To create maximin latin hypercube:
n<-16 # set length 16
library(lhs) #(This is an old package and has been known to crash)
# set empty vectors and distances
maxmind5<-0
maxmind4<-0
maxmind3<-0
maxmind2<-0
maxmind1<-0
superr5<-c()
superr4<-c()
superr3<-c()
superr2<-c()
superr1<-c()

while(maxmind1<0.25){ # set the minimum value you wish to have between each point (for plot took 0.20, but this takes some time)
  
  rl <- randomLHS(n,2) #create a random LHS
  
  d=min(dist(rl)) # find the minimum distance between all points
  # Then compare this minimum distance to the five previously saved. If it is larger than distance maxmind1, this LHS will become superr1, and push all
  # previous LHSs and their distances back one place, saving the 5 vectors with the largest minimum distances between points.
  if(d>maxmind5){
    if(d>maxmind4){
      if(d>maxmind3){
        if(d>maxmind2){
          if(d>maxmind1){
            superr5<-superr4
            superr4<-superr3
            superr3<-superr2
            superr2<-superr1
            
            maxmind5<-maxmind4
            maxmind4<-maxmind3
            maxmind3<-maxmind2
            maxmind2<-maxmind1
            
            maxmind1<-d
            superr1<-rl
          }
          
          else{
            superr5<-superr4
            superr4<-superr3
            superr3<-superr2
            
            maxmind5<-maxmind4
            maxmind4<-maxmind3
            maxmind3<-maxmind2
            
            maxmind2<-d
            superr2<-rl
          }
        }
        
        else{
          superr5<-superr4
          superr4<-superr3
          
          maxmind5<-maxmind4
          maxmind4<-maxmind3
          
          maxmind3<-d
          superr3<-rl
        }
      }
      else{
        superr5<-superr4
        maxmind5<-maxmind4
        
        maxmind4<-d
        superr4<-rl
      }
    }
    else{
      maxmind5<-d
      superr5<-rl
    }
  }
}
# This was then altered to swap min for max and reversed inequalities to find a minimax "bad" LHS.
# Both these were then simply plotted first column against the second.

# Plotting the maximin
plot(superr2, lwd=6,xlim = c(0,1),ylim=c(0,1),xlab = expression('x'[1]), ylab = expression('x'[2]))

##################
# Figure 3.5
##################
x <- seq(1.2,3.4,0.01)
y <- seq(1.2,3.4,0.01)
xy <- expand.grid(x,y)
f <- function(x){cos(0.5*x[,1])+sin(x[,1]*x[,2])}
B0 = 0.25
sigu = 1.5
theta = 0.8
csjs<-cbind(superr1[,1]*2.2+1.2,superr1[,2]*2.2+1.2)
js<-csjs # obviously due to its random nature the superr1 we used will be different to any other.
fp2(actual=TRUE,expectation=TRUE,SD=TRUE,implausibility=FALSE,SK=TRUE,Levelss=seq(-1.2,2,0.2),ssave=ss,nahme = 'cossinmmlh')



##################
# Figure 3.6
##################
x <- seq(0,2,0.025)
y <- seq(-1,1,0.025)
xy <- expand.grid(x,y)
f <- function(x){x[,2]*sin(4*x[,1])}
xjs <- seq(0.0,1.8,1.8/5)
yj1 <- rep(-1,6)
yj2 <- rep(0,6)
yj3 <- rep(1,6)
yjs <- c(yj1,yj2,yj3)
js <- cbind(xjs,yjs)
B0 = 0
sigu = 1.5
theta = 0.6
fp2(actual=FALSE,expectation=TRUE,SD=FALSE,implausibility=FALSE,SK=FALSE,Levelss=seq(-1.1,1.1,0.1),ssave=ss,nahme = 'goodysing')



##################
# Figure 3.7
##################
# We manipulate the function for a maximin LHS to achieve a min average sd LHS
# Again set empty vectors and sds
sumv5<-0
sumv4<-0
sumv3<-0
sumv2<-0
sumv1<-0
superjs5<-c()
superjs4<-c()
superjs3<-c()
superjs2<-c()
superjs1<-c()

n=16
xy2<-cbind(xy[,1],xy[,2])
for(j in seq(0,100,1)){ # takes 100 runs (we took 1 million ~ 1hr)
  
  rl <- randomLHS(n,2) # Take a random LHS
  js <- cbind(2.2*rl[,1]+1.2,2.2*rl[,2]+1.2) # Transform it to the plane
  
  covfD <- function(x){
    varf*exp(-(as.matrix(dist(rbind(x,js))))[(1+dim(x)[1]):(dim(x)[1]+dim(js)[1]),1:dim(x)[1]]^2/theta^2)
  }
  varD <- covfD(js)
  invvarD <- solve(varD)
  
  varDf <- function(x){
    A<-covfD(x)
    B<-invvarD%*%A
    varf-colSums(A * B)
  }
  z<-c()
  for (i in 0:floor((dim(xy2)[1]/100)-1)) {
    z<-c(z,c(varDf(xy2[(i*100):(i*100+99),])))
  }
  q<-dim(xy2)[1]-100*floor((dim(xy2)[1]/100)-1)-100
  z<-c(z,c(varDf(xy2[(dim(xy2)[1]-q):(dim(xy2)[1]),])))
  v=sum(z^0.5) # we then compare its total standard deviation (no different from comparing means)
  if(v<sumv5){
    if(v<sumv4){
      if(v<sumv3){
        if(v<sumv2){
          if(v<sumv1){
            superjs5<-superjs4
            superjs4<-superjs3
            superjs3<-superjs2
            superjs2<-superjs1
            
            sumv5<-sumv4
            sumv4<-sumv3
            sumv3<-sumv2
            sumv2<-sumv1
            
            sumv1<-v
            superjs1<-js
          }
          
          else{
            superjs5<-superjs4
            superjs4<-superjs3
            superjs3<-superjs2
            
            sumv5<-sumv4
            sumv4<-sumv3
            sumv3<-sumv2
            
            sumv2<-v
            superjs2<-js
          }
        }
        
        else{
          superjs5<-superjs4
          superjs4<-superjs3
          
          sumv5<-sumv4
          sumv4<-sumv3
          
          sumv3<-v
          superjs3<-js
        }
      }
      else{
        superjs5<-superjs4
        sumv5<-sumv4
        
        sumv4<-v
        superjs4<-js
      }
    }
    else{
      sumv5<-v
      superjs5<-js
    }
  }
}
# Leaving us with superjs1 as the minimum mean sd LHS, which we now plot

x <- seq(1.2,3.4,0.01)
y <- seq(1.2,3.4,0.01)
xy <- expand.grid(x,y)
f <- function(x){cos(0.5*x[,1])+sin(x[,1]*x[,2])}
B0 = 0.25
sigu = 1.5
theta = 0.8
js<-superjs1 # again different due to its random nature.
fp2(actual=TRUE,expectation=TRUE,SD=TRUE,implausibility=FALSE,SK=TRUE,Levelss=seq(-1.2,2,0.2),ssave=ss,nahme = 'cossinmmlh')



############
# Figure 3.8
############
# First we define a function to calculate average sd
avsd<-function(js){
  varf <- sigu^2
  xy2<-cbind(xy[,1],xy[,2])
  covfD <- function(x){
    varf*exp(-(as.matrix(dist(rbind(x,js))))[(1+dim(x)[1]):(dim(x)[1]+dim(js)[1]),1:dim(x)[1]]^2/theta^2)
  }
  varD <- covfD(js)
  invvarD <- solve(varD)
  
  varDf <- function(x){
    A<-covfD(x)
    B<-invvarD%*%A
    varf-colSums(A * B)
  }
  z<-c()
  for (i in 0:floor((dim(xy2)[1]/100)-1)) {
    z<-c(z,c(varDf(xy2[(i*100):(i*100+99),])^0.5))
  }
  q<-dim(xy2)[1]-100*floor((dim(xy2)[1]/100)-1)-100
  z<-c(z,c(varDf(xy2[(dim(xy2)[1]-q):(dim(xy2)[1]),])^0.5))
  z[which(is.na(z)==TRUE)]<-0 # for some designs single values came up with errors so we set them to 0. As this was only 1 value we believe
                              # that the effect is insignificant
  v=sum(z) # sum total sd
  return(v/nrow(xy2)) # divide by number of coordinates to get mean sd per coordinate
}

# Then create a function that takes a vector not matrix as input to calculate avsd(js)
sdtest<-function(c){
  dim(c)<-c(16,2)
  js<-c
  avsd(js)
}
c<-c(superjs1) #starting with previous design, set it as a vector
sdtest(c) # test its sd
rop<-optim(c,sdtest) # optimise to min sd
ct<-rop$par # sets ct as this optimised c
dim(ct)<-c(16,2) # set the dimensions back to normal
sdtest(ct) # see optimisation has worked

# Then repeat the process 3 times (we did this 10000)
for(i in 1:3){
  rop<-optim(ct,sdtest)
  ct<-rop$par
}
sdtest(ct) # see improvement
dim(ct)<-c(16,2)
optsjs<-ct # store this optimised js, and then plot

x <- seq(1.2,3.4,0.01)
y <- seq(1.2,3.4,0.01)
xy <- expand.grid(x,y)
f <- function(x){cos(0.5*x[,1])+sin(x[,1]*x[,2])}
B0 = 0.25
sigu = 1.5
theta = 0.8
js<-optsjs
fp2(actual=TRUE,expectation=TRUE,SD=TRUE,implausibility=FALSE,SK=TRUE,Levelss=seq(-1.2,2,0.2),ssave=ss,nahme = 'cossinopt')



############
# Figure 3.9
############
# First we need a function that finds a new acceptable point to simulate
# fnF takes the original js, closeness: the minimum distance between new points, and ILEVEL: the maximum implausibility value of a new point.
fnF<-function(js,closeness,ILEVEL){
  # Again emulate as normal
  D <- f(js)
  Ef <- B0
  varf <- sigu^2
  ED <- rep(Ef,length(D))
  covfD <- function(x){
    varf*exp(-(as.matrix(dist(rbind(x,js))))[1,-1]^2/theta^2)
  }
  varD <- c()
  for (i in seq(1,dim(js)[1],1)) {
    varD <- c(varD,covfD(js[i,]))
  }
  dim(varD) <- c(dim(js)[1],dim(js)[1])
  invvarD <- solve(varD)
  EDf <- function(x){Ef+covfD(x)%*%invvarD%*%(D-ED)}
  varDf <- function(x){varf-covfD(x)%*%invvarD%*%c(covfD(x))}
  # Define implausibility^2 
  I2<- function(x){((EDf(x)-valuu)^2)/(varDf(x)+varee)}
  
  g<-range(x)
  m<-0 # This indicates to stop when we are succesful
  l<-length(js[,1])
  while(m<1){
    k<-runif(2,g[1],g[2]) # First take a random point on the range of this plane
    kc<-0 # This is a 'score' for if the point is too close to any point in js
    if(I2(k)<ILEVEL^2){ # If k's implausibility is too high we reject it
      trf<-(as.matrix(dist(rbind(k,js)))[1,-1]<closeness) # Set a matrix that measures the distance between x and each point, returning true if it is too close.
      for(i in trf){
        if(i==FALSE){kc<-kc+1} # For every point that isn't too close, as 1 point to kc
      }
    }
    if(kc==l){ # If k was far enough from every point, we have kc equal to the number of points, if this is the case we accept k
      m<-m+1 # Terminate the loop when successful
    }
  }
  k # Return k
}

# This function allows us to add one good point at a time, so we now define a function to use this to add n points at a time
newjs<-function(n,closeness,ILEVEL=4){
  newjs2<-js
  m<-0
  while(m<n){
    ks<-fnF(newjs2,closeness,ILEVEL) # Add a new point
    newjs2<-rbind(newjs2,ks) # Redefine the input
    m<-m+1 
  }
  newjs2 # Return new set of js plus n new points checked against js and each other
}

# Now to make the plots
x <- seq(1.2,3.4,0.01)
y <- seq(1.2,3.4,0.01)
xy <- expand.grid(x,y)
f <- function(x){cos(0.5*x[,1])+sin(x[,1]*x[,2])}
B0 = 0.25
sigu = 1.5
theta = 0.8
js<-optsjs
varee=(0.01*3.2)^2 # Set Var[epsilon]=0.001024
valuu=-0.5 # Set the value we try to match as -0.5
fp2(actual=TRUE,expectation=FALSE,SD=FALSE,implausibility=TRUE,SK=FALSE,Levelss=seq(-1.2,2,0.2),ssave=ss,nahme = 'cossinI')
js<-newjs(3,0.3,3) # Find three new points with implausibility < 3 and closeness > 0.3
fp2(actual=FALSE,expectation=FALSE,SD=FALSE,implausibility=TRUE,SK=FALSE,Levelss=seq(-1.2,2,0.2),ssave=ss,nahme = 'cossinI2')
js<-newjs(3,0.3,3)
fp2(actual=FALSE,expectation=FALSE,SD=FALSE,implausibility=TRUE,SK=FALSE,Levelss=seq(-1.2,2,0.2),ssave=ss,nahme = 'cossinI2')
js<-newjs(3,0.1,3) # We allowed closer values for the last plot as there was a lot less space to choose from
fp2(actual=FALSE,expectation=FALSE,SD=FALSE,implausibility=TRUE,SK=FALSE,Levelss=seq(-1.2,2,0.2),ssave=ss,nahme = 'cossinI3')



##############
# Calculations
##############
# We used function avsd(js) to calculate the mean standard deviation for each design js.
# We define function SK(js) to calculate the mean diagnostic value for each design js.
SK<-function(js){
  xy2<-cbind(xy[,1],xy[,2])
  D <- f(js)
  Ef <- B0
  varf <- sigu^2
  ED <- rep(Ef,length(D))
  covfD <- function(x){
    varf*exp(-(as.matrix(dist(rbind(x,js))))[(1+dim(x)[1]):(dim(x)[1]+dim(js)[1]),1:dim(x)[1]]^2/theta^2)
  }
  varD <- covfD(js)
  invvarD <- solve(varD)
  EDf <- function(x){
    A<-covfD(x)
    C<-rep(D-ED,dim(x)[1])
    dim(C)<-c(length(D-ED),dim(x)[1])
    B<-invvarD%*%C
    Ef+colSums(A * B)
  }
  varDf <- function(x){
    A<-covfD(x)
    B<-invvarD%*%A
    varf-colSums(A * B)
  }
  SC<- function(x){(EDf(x)-f(x))/((varDf(x))^0.5+0.000000000001)}
  # To calculate average (EDf(x)-f(x)) value instead, we uncomment the line below
  #SC<- function(x){EDf(x)-f(x)}
  z<-c()
  for (i in 0:floor((dim(xy2)[1]/100)-1)) {
    z<-c(z,c(SC(xy2[(i*100):(i*100+99),])))
  }
  q<-dim(xy2)[1]-100*floor((dim(xy2)[1]/100)-1)-100
  z<-c(z,c(SC(xy2[(dim(xy2)[1]-q):(dim(xy2)[1]),])))
  z[which(is.na(z)==TRUE)]<-0 # for some designs single values came up with errors so we set them to 0. As this was only 1 value we believe
  # that the effect is insignificant
  sum(abs(z))/length(z)
}
###
x <- seq(1.2,3.4,0.01)
y <- seq(1.2,3.4,0.01)
xy <- expand.grid(x,y)
f <- function(x){cos(0.5*x[,1])+sin(x[,1]*x[,2])}
xjs <- seq(1.2,3.4,2.2/3)
yjs <- seq(1.2,3.4,2.2/3)
js <- as.matrix(expand.grid(xjs,yjs))
B0 = 0.25
sigu = 1.5
theta = 0.8

SK(js)
#0.383904
avsd(js)
#0.3503486

js<-csjs
SK(js)
#0.3522398
avsd(js)
#0.2326086

js<-superjs1
SK(js)
#0.5544466
avsd(js)
#0.2265681

js<-optsjs
SK(js)
#0.2402472
avsd(js)
#0.2077678





####################################################################
#
#       Chapter 4
#
####################################################################
# For this chapter we define three new functions, similar to fp2, but now taking into account the new
# equations for E, Cov and Var according to the different boundary conditions.
# First fpB2 emulates with a single known boundary along the left y-axis (This can be easily moved in the code) and no design
fpB2<-function(actualplot,expectation,SD,SK,Levelss,ssave,nahme){
  xy2<-cbind(xy[,1],xy[,2])
  if(actualplot==TRUE){
    z<-f(xy2)
    dim(z) <- c(length(x),length(y))
    
    if(ssave==TRUE){pdf(file=file.path("C:","Users","tcvin","OneDrive","Work","DIS","pics","core",paste(nahme,"a",".pdf",sep="")))}
    plot.new()
    plot.window(xlim = c(1.2,3.4),ylim = c(1.2,3.4),asp=NA)
    .filled.contour(x,y, z,levels = Levelss, col = rainbow(16))
    if(ssave==TRUE){dev.off()}
    
  }
  K<-x[1] # Here we define K to run along the first x_1 value x[1], which we can easily change
  Ef <- B0
  varf <- sigu^2
  EKf <- function(x){ # E_K[f(x)]
    a<-x[,1]-K # Define a as the distance from x to K
    Ef+exp(-abs(a)^2/theta^2)*(f(cbind(K,x[,2]))-Ef)
  }
  varKf <- function(x){ # Var_K[f(x)]
    a<-x[,1]-K
    varf*(1-exp(-2*abs(a)^2/theta^2))
  }
  if(expectation==TRUE){
    z<-EKf(xy2) # As no dist function is used, emulation of very large matrices is still efficient
    dim(z) <- c(length(x),length(y))
    if(ssave==TRUE){pdf(file=file.path("C:","Users","tcvin","OneDrive","Work","DIS","pics","core",paste(nahme,"e",".pdf",sep="")))}
    plot.new()
    plot.window(xlim = c(1.2,3.4),ylim = c(1.2,3.4),asp=NA)
    .filled.contour(x,y, z,levels = Levelss, col = rainbow(16))
    if(ssave==TRUE){dev.off()}
  }
  if(SD==TRUE){
    z<-varKf(xy2)^0.5
    dim(z) <- c(length(x),length(y))
    if(ssave==TRUE){pdf(file=file.path("C:","Users","tcvin","OneDrive","Work","DIS","pics","core",paste(nahme,"sd",".pdf",sep="")))}
    plot.new()
    plot.window(xlim = c(1.2,3.4),ylim = c(1.2,3.4),asp=NA)
    .filled.contour(x,y, z,levels = seq(0,1.5,0.1), col = terrain.colors(15))
    if(ssave==TRUE){dev.off()}
  }
  if(SK==TRUE){
    SC<- function(x){(EKf(x)-f(x))/((varKf(x))^0.5+0.000000000001)}
    z<-SC(xy2)
    dim(z) <- c(length(x),length(y))
    if(ssave==TRUE){pdf(file=file.path("C:","Users","tcvin","OneDrive","Work","DIS","pics","core",paste(nahme,"sk",".pdf",sep="")))}
    plot.new()
    plot.window(xlim = c(1.2,3.4),ylim = c(1.2,3.4),asp=NA)
    .filled.contour(x,y, z,levels = seq(-2,2,0.1), col = colorRampPalette(c('red','chartreuse','navy'),1)(41))
    if(ssave==TRUE){dev.off()}
  }
}

# fpBD2 now emulates with single K and some design D defined by js
fpBD2<-function(actualplot,expectation,SD,SK,Levelss,ssave,nahme){
  xy2<-cbind(xy[,1],xy[,2])
  if(actualplot==TRUE){
    z<-f(xy2)
    dim(z) <- c(length(x),length(y))
    
    if(ssave==TRUE){pdf(file=file.path("C:","Users","tcvin","OneDrive","Work","DIS","pics","core",paste(nahme,"a",".pdf",sep="")))}
    plot.new()
    plot.window(xlim = c(1.2,3.4),ylim = c(1.2,3.4),asp=NA)
    .filled.contour(x,y, z,levels = Levelss, col = rainbow(16))
    points(js[,1],js[,2],lwd=6)
    if(ssave==TRUE){dev.off()}
    
  }
  K <- x[1]
  Ef <- B0
  varf <- sigu^2
  EKf <- function(x){
    a<-x[,1]-K
    Ef+exp(-abs(a)^2/theta^2)*(f(cbind(K,x[,2]))-Ef)
  }
  varKf <- function(x){
    a<-x[,1]-K
    varf*(1-exp(-2*abs(a)^2/theta^2))
  }
  R<-function(a1,a2){ # Define R as required
    exp(-abs(a1-a2)^2/theta^2)-exp((-abs(a1)^2-abs(a2)^2)/theta^2)
  }
  covKfD <- function(x1){ # We define Cov_K[f(x),D], where D is predefined, and x (here being x1) is the only input
    a1<-x1[,1]-K
    a2<-js[,1]-K
    # We calculate q = -abs(xK2-x'K2)^2/theta^2, where input x1=x and design js=x'
    # Say we input n coordinates as x1, and compare to 16 points js as in the report. The first matrix repeats xK2 16 times for each of the n coordinates.
    # The second matrix repeats the 16 x'K2 values for each point in js, and does this for n rows. We then simply subtract these two matrices, take the
    # absolute, square and divide by theta.
    q<- -abs((matrix(data = rep(x1[,2],length(js[,2])), nrow = length(x1[,2])))-t(matrix(data = rep(js[,2],length(x1[,2])), ncol = length(x1[,2]))))^2/theta^2
    #varf*R(a1,a2)*exp(-abs(x1[,2]-x2)^2/theta^2)
    eq<-exp(q) # exp[-abs(xK2-x'K2)^2/theta^2]
    aa<-expand.grid(a2,a1) # a and a' value for each pair of x coordinate and point in js
    RR<-R(aa[,2],aa[,1]) # R(a,a') value for each pair of x coordinate and point in js
    dim(RR)<-c(length(js[,2]),length(x1[,2])) # Make right dimensions then calculate full equation.
    RR<-t(RR)
    varf*RR*eq
  }
  D <- f(js)
  varKD <- covKfD(js)
  invvarKD <- solve(varKD)
  EDKf <- function(x){
    A<-t(covKfD(x))
    C<-(matrix(data=rep(D-EKf(js),dim(x)[1]),ncol=dim(x)[1]))
    B<-invvarKD%*%C
    EKf(x)+colSums(A * B)
  }
  varDKf <- function(x){
    A<-t(covKfD(x))
    B<-invvarKD%*%A
    varKf(x)-colSums(A * B)
  }
  if(expectation==TRUE){
    z<-EDKf(xy2)
    z[which(z< -1.2)]<- 1*-1.19
    dim(z) <- c(length(x),length(y))
    if(ssave==TRUE){pdf(file=file.path("C:","Users","tcvin","OneDrive","Work","DIS","pics","core",paste(nahme,"e",".pdf",sep="")))}
    plot.new()
    plot.window(xlim = c(1.2,3.4),ylim = c(1.2,3.4),asp=NA)
    .filled.contour(x,y, z,levels = Levelss, col = rainbow(16))
    points(js[,1],js[,2],lwd=6)
    if(ssave==TRUE){dev.off()}
  }
  if(SD==TRUE){
    z<-varDKf(xy2)^0.5
    dim(z) <- c(length(x),length(y))
    if(ssave==TRUE){pdf(file=file.path("C:","Users","tcvin","OneDrive","Work","DIS","pics","core",paste(nahme,"sd",".pdf",sep="")))}
    plot.new()
    plot.window(xlim = c(1.2,3.4),ylim = c(1.2,3.4),asp=NA)
    .filled.contour(x,y, z,levels = seq(0,1,0.05), col = terrain.colors(15))
    points(js[,1],js[,2],lwd=6)
    if(ssave==TRUE){dev.off()}
  }
  if(SK==TRUE){
    SC<- function(x){(EDKf(x)-f(x))/((varDKf(x))^0.5+0.000000000001)}
    z<-SC(xy2)
    dim(z) <- c(length(x),length(y))
    if(ssave==TRUE){pdf(file=file.path("C:","Users","tcvin","OneDrive","Work","DIS","pics","core",paste(nahme,"sk",".pdf",sep="")))}
    
    plot.new()
    plot.window(xlim = c(1.2,3.4),ylim = c(1.2,3.4),asp=NA)
    .filled.contour(x,y, z,levels = seq(-2,2,0.1), col = colorRampPalette(c('red','chartreuse','navy'),1)(41))
    points(js[,1],js[,2],lwd=6)
    if(ssave==TRUE){dev.off()}
  }
}

# fpBPerp2 emulates with perpendicular boundaries, no design
fpBPerp2<-function(actualplot,expectation,SD,SK,Levelss,ssave,nahme){
  xy2<-cbind(xy[,1],xy[,2])
  w<-xy2[1,1] # Here w marks the point L and K meet, thus defining both. Chosen as along the left and bottom axes.
  if(actualplot==TRUE){
    z<-f(xy2)
    dim(z) <- c(length(x),length(y))
    
    if(ssave==TRUE){pdf(file=file.path("C:","Users","tcvin","OneDrive","Work","DIS","pics","core",paste(nahme,"e",".pdf",sep="")))}
    plot.new()
    plot.window(xlim = c(1.2,3.4),ylim = c(1.2,3.4),asp=NA)
    .filled.contour(x,y, z,levels = Levelss, col = rainbow(16))
    if(ssave==TRUE){dev.off()}
    
  }
  K <- f(cbind(rep(x[1],length(x)),y))
  Ef <- B0
  varf <- sigu^2
  r<-function(a){ # We define this to save time retyping
    exp(-abs(a)^2/theta^2)
  }
  nab<-function(x){ # nabla
    f(x)-Ef
  }
  ELKf <- function(x){ # E_LK[f(x)]
    a<-x[,1]-w
    b<-x[,2]-w
    xK<-cbind(w,x[,2])
    xL<-cbind(x[,1],w)
    xLK<-cbind(w,xL[,2])
    Ef+r(a)*nab(xK)+r(b)*nab(xL)-r(a)*r(b)*nab(xLK)
  }
  varLKf <- function(x){ # Var_LK[f(x)]
    a<-x[,1]-w
    b<-x[,2]-w
    varf*R(a,a)*R(b,b)
  }
  if(expectation==TRUE){
    z<-ELKf(xy2)
    dim(z) <- c(length(x),length(y))
    if(ssave==TRUE){pdf(file=file.path("C:","Users","tcvin","OneDrive","Work","DIS","pics","core",paste(nahme,"e",".pdf",sep="")))}
    plot.new()
    plot.window(xlim = c(1.2,3.4),ylim = c(1.2,3.4),asp=NA)
    .filled.contour(x,y, z,levels = Levelss, col = rainbow(16))
    if(ssave==TRUE){dev.off()}
  }
  if(SD==TRUE){
    z<-varLKf(xy2)^0.5
    dim(z) <- c(length(x),length(y))
    if(ssave==TRUE){pdf(file=file.path("C:","Users","tcvin","OneDrive","Work","DIS","pics","core",paste(nahme,"sd",".pdf",sep="")))}
    plot.new()
    plot.window(xlim = c(1.2,3.4),ylim = c(1.2,3.4),asp=NA)
    .filled.contour(x,y, z,levels = seq(0,1.5,0.1), col = terrain.colors(15))
    if(ssave==TRUE){dev.off()}
  }
  if(SK==TRUE){
    SC<- function(x){(ELKf(x)-f(x))/((varLKf(x))^0.5+0.000000000001)}
    z<-SC(xy2)
    dim(z) <- c(length(x),length(y))
    if(ssave==TRUE){pdf(file=file.path("C:","Users","tcvin","OneDrive","Work","DIS","pics","core",paste(nahme,"sk",".pdf",sep="")))}
    plot.new()
    plot.window(xlim = c(1.2,3.4),ylim = c(1.2,3.4),asp=NA)
    .filled.contour(x,y, z,levels = seq(-2,2,0.1), col = colorRampPalette(c('red','chartreuse','navy'),1)(41))
    if(ssave==TRUE){dev.off()}
  }
}

# fpBPerpD2 emulates with perpendicular boundaries and design js
fpBPerpD2<-function(actualplot,expectation,SD,SK,Levelss,ssave,nahme){
  xy2<-cbind(xy[,1],xy[,2])
  w<-xy2[1,1]
  if(actualplot==TRUE){
    z<-f(xy2)
    dim(z) <- c(length(x),length(y))
    
    if(ssave==TRUE){pdf(file=file.path("C:","Users","tcvin","OneDrive","Work","DIS","pics","core",paste(nahme,"e",".pdf",sep="")))}
    plot.new()
    plot.window(xlim = c(1.2,3.4),ylim = c(1.2,3.4),asp=NA)
    .filled.contour(x,y, z,levels = Levelss, col = rainbow(16))
    points(js[,1],js[,2],lwd=6)
    if(ssave==TRUE){dev.off()}
    
  }
  Ef <- B0
  varf <- sigu^2
  r<-function(a){
    exp(-abs(a)^2/theta^2)
  }
  nab<-function(x){
    f(x)-Ef
  }
  ELKf <- function(x){
    a<-x[,1]-w
    b<-x[,2]-w
    xK<-cbind(w,x[,2])
    xL<-cbind(x[,1],w)
    xLK<-cbind(w,xL[,2])
    Ef+r(a)*nab(xK)+r(b)*nab(xL)-r(a)*r(b)*nab(xLK)
  }
  R<-function(a1,a2){
    exp(-abs(a1-a2)^2/theta^2)-exp((-abs(a1)^2-abs(a2)^2)/theta^2)
  }
  varLKf <- function(x){
    a<-x[,1]-w
    b<-x[,2]-w
    varf*R(a,a)*R(b,b)
  }
  covLKfD <- function(x1){ # Cov_LK[f(x),D] similarly to before, only compares x with js
    a1<-x1[,1]-xy2[1,1] # a
    a2<-js[,1]-xy2[1,1] # a'
    b1<-x1[,2]-xy2[1,1] # b
    b2<-js[,2]-xy2[1,1] # b'
    aa<-expand.grid(a2,a1)
    bb<-expand.grid(b2,b1)
    RRa<-R(aa[,2],aa[,1])
    dim(RRa)<-c(length(js[,2]),length(x1[,2]))
    RRa<-t(RRa) # R(a,a')
    RRb<-R(bb[,2],bb[,1])
    dim(RRb)<-c(length(js[,2]),length(x1[,2]))
    RRb<-t(RRb) # R(b,b')
    varf*RRa*RRb
  }
  D <- f(js)
  varLKD <- covLKfD(js)
  invvarLKD <- solve(varLKD)
  EDLKf <- function(x){ # E_DLK[f(x)]
    A<-t(covLKfD(x))
    C<-(matrix(data=rep(D-ELKf(js),dim(x)[1]),ncol=dim(x)[1]))
    B<-invvarLKD%*%C
    ELKf(x)+colSums(A * B)
  }
  varDLKf <- function(x){ # Var_DLK[f(x)]
    A<-t(covLKfD(x))
    B<-invvarLKD%*%A
    varLKf(x)-colSums(A * B)
  }
  if(expectation==TRUE){
    z<-EDLKf(xy2)
    z[which(z< -1.2)]<- 1*-1.19
    dim(z) <- c(length(x),length(y))
    if(ssave==TRUE){pdf(file=file.path("C:","Users","tcvin","OneDrive","Work","DIS","pics","core",paste(nahme,"e",".pdf",sep="")))}
    plot.new()
    plot.window(xlim = c(1.2,3.4),ylim = c(1.2,3.4),asp=NA)
    .filled.contour(x,y, z,levels = Levelss, col = rainbow(16))
    points(js[,1],js[,2],lwd=6)
    if(ssave==TRUE){dev.off()}
  }
  if(SD==TRUE){
    z<-varDLKf(xy2)^0.5
    dim(z) <- c(length(x),length(y))
    if(ssave==TRUE){pdf(file=file.path("C:","Users","tcvin","OneDrive","Work","DIS","pics","core",paste(nahme,"sd",".pdf",sep="")))}
    plot.new()
    plot.window(xlim = c(1.2,3.4),ylim = c(1.2,3.4),asp=NA)
    .filled.contour(x,y, z,levels = seq(0,1,0.05), col = terrain.colors(20))
    points(js[,1],js[,2],lwd=6)
    if(ssave==TRUE){dev.off()}
  }
  if(SK==TRUE){
    SC<- function(x){(EDLKf(x)-f(x))/((varDLKf(x))^0.5+0.000000000001)}
    z<-SC(xy2)
    dim(z) <- c(length(x),length(y))
    if(ssave==TRUE){pdf(file=file.path("C:","Users","tcvin","OneDrive","Work","DIS","pics","core",paste(nahme,"sk",".pdf",sep="")))}
    plot.new()
    plot.window(xlim = c(1.2,3.4),ylim = c(1.2,3.4),asp=NA)
    .filled.contour(x,y, z,levels = seq(-2,2,0.1), col = colorRampPalette(c('red','chartreuse','navy'),1)(41))
    points(js[,1],js[,2],lwd=6)
    if(ssave==TRUE){dev.off()}
  }
}

# fpBPara2 emulates with perpendicular boundaries and no design
fpBPara2<-function(actualplot,expectation,SD,SK,Levelss,ssave,nahme){
  xy2<-cbind(xy[,1],xy[,2])
  K<-xy2[1,1] # Choses K and L as boundaries on the left and right
  L<-xy2[length(x),1]
  if(actualplot==TRUE){
    z<-f(xy2)
    dim(z) <- c(length(x),length(y))
    
    if(ssave==TRUE){pdf(file=file.path("C:","Users","tcvin","OneDrive","Work","DIS","pics","core",paste(nahme,"e",".pdf",sep="")))}
    plot.new()
    plot.window(xlim = c(1.2,3.4),ylim = c(1.2,3.4),asp=NA)
    .filled.contour(x,y, z,levels = Levelss, col = rainbow(16))
    if(ssave==TRUE){dev.off()}
    
  }
  Ef <- B0
  varf <- sigu^2
  r<-function(a){
    exp(-abs(a)^2/theta^2)
  }
  nab<-function(x){
    f(x)-Ef
  }
  R<-function(a1,a2){
    exp(-abs(a1-a2)^2/theta^2)-exp((-abs(a1)^2-abs(a2)^2)/theta^2)
  }
  RrR<-function(a,b,c){
    (r(a)-r(b)*r(c))/(1-r(c)^2)
  }
  ELKf <- function(x){ # E_LK[f(x)]
    a<-x[,1]-K
    b<-L-x[,1]
    c<-L-K
    xK<-cbind(K,x[,2])
    xL<-cbind(L,x[,2])
    Ef+RrR(a,b,c)*nab(xK)+RrR(b,a,c)*nab(xL)
  }
  varLKf <- function(x){ # Var_LK[f(x)]
    a<-x[,1]-K
    b<-L-x[,1]
    c<-L-K
    (varf/(1-r(c)^2))*(1-r(c)^2-r(a)^2-r(b)^2+2*r(a)*r(b)*r(c))
  }
  if(expectation==TRUE){
    z<-ELKf(xy2)
    dim(z) <- c(length(x),length(y))
    if(ssave==TRUE){pdf(file=file.path("C:","Users","tcvin","OneDrive","Work","DIS","pics","core",paste(nahme,"e",".pdf",sep="")))}
    plot.new()
    plot.window(xlim = c(1.2,3.4),ylim = c(1.2,3.4),asp=NA)
    .filled.contour(x,y, z,levels = Levelss, col = rainbow(16))
    if(ssave==TRUE){dev.off()}
  }
  if(SD==TRUE){
    z<-varLKf(xy2)^0.5
    dim(z) <- c(length(x),length(y))
    if(ssave==TRUE){pdf(file=file.path("C:","Users","tcvin","OneDrive","Work","DIS","pics","core",paste(nahme,"sd",".pdf",sep="")))}
    plot.new()
    plot.window(xlim = c(1.2,3.4),ylim = c(1.2,3.4),asp=NA)
    .filled.contour(x,y, z,levels = seq(0,1.5,0.1), col = terrain.colors(15))
    if(ssave==TRUE){dev.off()}
  }
  if(SK==TRUE){
    SC<- function(x){(ELKf(x)-f(x))/((varLKf(x))^0.5+0.000000000001)}
    z<-SC(xy2)
    dim(z) <- c(length(x),length(y))
    if(ssave==TRUE){pdf(file=file.path("C:","Users","tcvin","OneDrive","Work","DIS","pics","core",paste(nahme,"sk",".pdf",sep="")))}
    plot.new()
    plot.window(xlim = c(1.2,3.4),ylim = c(1.2,3.4),asp=NA)
    .filled.contour(x,y, z,levels = seq(-2,2,0.1), col = colorRampPalette(c('red','chartreuse','navy'),1)(41))
    if(ssave==TRUE){dev.off()}
  }
}

# fpBParaD2 emulates with perpendicular boundaries and design js
fpBParaD2<-function(actualplot,expectation,SD,SK,Levelss,ssave,nahme){
  xy2<-cbind(xy[,1],xy[,2])
  K<-xy2[1,1]
  L<-xy2[length(x),1]
  if(actualplot==TRUE){
    z<-f(xy2)
    dim(z) <- c(length(x),length(y))
    
    if(ssave==TRUE){pdf(file=file.path("C:","Users","tcvin","OneDrive","Work","DIS","pics","core",paste(nahme,"e",".pdf",sep="")))}
    plot.new()
    plot.window(xlim = c(1.2,3.4),ylim = c(1.2,3.4),asp=NA)
    .filled.contour(x,y, z,levels = Levelss, col = rainbow(16))
    points(js[,1],js[,2],lwd=6)
    if(ssave==TRUE){dev.off()}
    
  }
  Ef <- B0
  varf <- sigu^2
  r<-function(a){
    exp(-abs(a)^2/theta^2)
  }
  nab<-function(x){
    f(x)-Ef
  }
  R<-function(a1,a2){
    exp(-abs(a1-a2)^2/theta^2)-exp((-abs(a1)^2-abs(a2)^2)/theta^2)
  }
  RrR<-function(a,b,c){
    (r(a)-r(b)*r(c))/(1-r(c)^2)
  }
  ELKf <- function(x){
    a<-x[,1]-K
    b<-L-x[,1]
    c<-L-K
    xK<-cbind(K,x[,2])
    xL<-cbind(L,x[,2])
    Ef+RrR(a,b,c)*nab(xK)+RrR(b,a,c)*nab(xL)
  }
  varLKf <- function(x){
    a<-x[,1]-K
    b<-L-x[,1]
    c<-L-K
    (varf/(1-r(c)^2))*(1-r(c)^2-r(a)^2-r(b)^2+2*r(a)*r(b)*r(c))
  }
  covLKfD <- function(x1){ # Cov_LK[f(x)]
    a1<-x1[,1]-K
    a2<-js[,1]-K
    c<-L-K
    aa<-expand.grid(a2,a1)
    RRa<-R(aa[,2],aa[,1])
    dim(RRa)<-c(length(js[,2]),length(x1[,2]))
    RRa<-t(RRa)
    ac<-cbind(aa[,2],c)
    RRac<-R(ac[,1],ac[,2])
    dim(RRac)<-c(length(js[,2]),length(x1[,2]))
    RRac<-t(RRac)
    RRca<-R(aa[,1],ac[,2])
    dim(RRca)<-c(length(js[,2]),length(x1[,2]))
    RRca<-t(RRca)
    q<- -abs(t(matrix(data = rep(js[,2],length(x1[,2])), ncol = length(x1[,2])))-(matrix(data = rep(x1[,2],length(js[,2])), nrow = length(x1[,2]))))^2/theta^2
    eq<-exp(q)
    varf*eq*(RRa-(RRac*RRca/R(c,c)))
  }
  D <- f(js)
  varLKD <- covLKfD(js)
  invvarLKD <- solve(varLKD)
  EDLKf <- function(x){ # E_DLK[f(x)]
    A<-t(covLKfD(x))
    C<-(matrix(data=rep(D-ELKf(js),dim(x)[1]),ncol=dim(x)[1]))
    B<-invvarLKD%*%C
    ELKf(x)+colSums(A * B)
  }
  varDLKf <- function(x){ # Var_DLK[f(x)]
    A<-t(covLKfD(x))
    B<-invvarLKD%*%A
    varLKf(x)-colSums(A * B)
  }
  if(expectation==TRUE){
    z<-EDLKf(xy2)
    dim(z) <- c(length(x),length(y))
    if(ssave==TRUE){pdf(file=file.path("C:","Users","tcvin","OneDrive","Work","DIS","pics","core",paste(nahme,"e",".pdf",sep="")))}
    plot.new()
    plot.window(xlim = c(1.2,3.4),ylim = c(1.2,3.4),asp=NA)
    .filled.contour(x,y, z,levels = Levelss, col = rainbow(16))
    points(js[,1],js[,2],lwd=6)
    if(ssave==TRUE){dev.off()}
  }
  if(SD==TRUE){
    z<-varDLKf(xy2)^0.5
    dim(z) <- c(length(x),length(y))
    if(ssave==TRUE){pdf(file=file.path("C:","Users","tcvin","OneDrive","Work","DIS","pics","core",paste(nahme,"sd",".pdf",sep="")))}
    plot.new()
    plot.window(xlim = c(1.2,3.4),ylim = c(1.2,3.4),asp=NA)
    .filled.contour(x,y, z,levels = seq(0,1,0.05), col = terrain.colors(20))
    points(js[,1],js[,2],lwd=6)
    if(ssave==TRUE){dev.off()}
  }
  if(SK==TRUE){
    SC<- function(x){(EDLKf(x)-f(x))/((varDLKf(x))^0.5+0.000000000001)}
    z<-SC(xy2)
    dim(z) <- c(length(x),length(y))
    if(ssave==TRUE){pdf(file=file.path("C:","Users","tcvin","OneDrive","Work","DIS","pics","core",paste(nahme,"sk",".pdf",sep="")))}
    plot.new()
    plot.window(xlim = c(1.2,3.4),ylim = c(1.2,3.4),asp=NA)
    .filled.contour(x,y, z,levels = seq(-2,2,0.1), col = colorRampPalette(c('red','chartreuse','navy'),1)(41))
    points(js[,1],js[,2],lwd=6)
    if(ssave==TRUE){dev.off()}
  }
}

#############
# FIGURES
#############
# In order to create the optimised designs, we first manipulated the code from figure 3.7, replacing the no boundary calculations of
# variance with the relevent known boundary calculations from above. Then we used the same code as in 3.8 to optimise, again altering
# sdtest with relevant variance caulations. To evidence how the plots worked we use the previously optimised js.

x <- seq(1.2,3.4,0.01)
y <- seq(1.2,3.4,0.01)
xy <- expand.grid(x,y)
f <- function(x){cos(0.5*x[,1])+sin(x[,1]*x[,2])}
B0 = 0.25
sigu = 1.5
theta = 0.8
js<-optsjs

#4.2
fpB2(actual=TRUE,expectation=TRUE,SD=TRUE,SK=TRUE,Levelss=seq(-1.2,2,0.2),ssave=ss,nahme = 'csB')
#4.3
fpBD2(actual=FALSE,expectation=TRUE,SD=TRUE,SK=TRUE,Levelss=seq(-1.2,2,0.2),ssave=ss,nahme = 'csBD')
#4.5
fpBPerp2(actual=FALSE,expectation=TRUE,SD=TRUE,SK=TRUE,Levelss=seq(-1.2,2,0.2),ssave=ss,nahme = 'csBPerp')
fpBPerpD2(actual=FALSE,expectation=TRUE,SD=TRUE,SK=TRUE,Levelss=seq(-1.2,2,0.2),ssave=ss,nahme = 'csBPerpD')
#4.6
fpBPara2(actual=FALSE,expectation=TRUE,SD=TRUE,SK=TRUE,Levelss=seq(-1.2,2,0.2),ssave=ss,nahme = 'csBPara')
fpBParaD2(actual=FALSE,expectation=TRUE,SD=TRUE,SK=TRUE,Levelss=seq(-1.2,2,0.2),ssave=ss,nahme = 'csBParaD')





####################################################################
#
#       Chapter 5
#
####################################################################
# Here we define the functions that take 'slices' of the 3D input space for a fixed x1 or x2 value.
# We assume for all plots that each input runs between 0 and 1, and that L and K run as in the report.

#TD1 takes a slice for fixed x1 value (input x1value) and emulates without any design
TD1<-function(actualplot,expectation,SD,SK,Lev,ssave,nahme,x1value){
  x <- seq(0,1,0.01)
  y <- seq(0,1,0.01)
  xy <- expand.grid(x,y)
  xy2<-cbind(x1value,xy[,1],xy[,2])
  w<-0 # Again define K and L to meet at (0,0)
  if(actualplot==TRUE){
    z<-f(xy2)
    dim(z) <- c(length(x),length(y))
    
    if(ssave==TRUE){png(file=file.path("C:","Users","tcvin","OneDrive","Work","DIS","pics","core",paste(nahme,"a",".png",sep="")))}
    plot.new()
    plot.window(xlim = c(0,1),ylim = c(0,1),asp=NA)
    .filled.contour(x,y, z,levels = Lev, col = rainbow(25))
    if(ssave==TRUE){dev.off()}
    
  }
  Ef <- B0
  varf <- sigu^2
  r<-function(a){
    exp(-abs(a)^2/theta^2)
  }
  R<-function(a1,a2){
    exp(-abs(a1-a2)^2/theta^2)-exp((-abs(a1)^2-abs(a2)^2)/theta^2)
  }
  nab<-function(x){
    f(x)-Ef
  }
  ELKf <- function(x){
    a<-x[,1]-w
    b<-x[,2]-w
    c<-x[,3]-w
    xK<-cbind(w,x[,2],x[,3])
    xL<-cbind(x[,1],w,w)
    xLK<-cbind(w,xL[,2],w)
    Ef+r(a)*nab(xK)+r(b)*r(c)*nab(xL)-r(a)*r(b)*r(c)*nab(xLK)
  }
  varLKf <- function(x){
    a<-x[,1]-w
    b<-x[,2]-w
    c<-x[,3]-w
    varf*(1-r(a)^2)*(1-r(b)^2*r(c)^2)
  }
  if(expectation==TRUE){
    z<-ELKf(xy2)
    dim(z) <- c(length(x),length(y))
    if(ssave==TRUE){png(file=file.path("C:","Users","tcvin","OneDrive","Work","DIS","pics","core",paste(nahme,"e",".png",sep="")))}
    plot.new()
    plot.window(xlim = c(0,1),ylim = c(0,1),asp=NA)
    .filled.contour(x,y, z,levels = Lev, col = rainbow(25))
    if(ssave==TRUE){dev.off()}
  }
  if(SD==TRUE){
    z<-varLKf(xy2)^0.5
    dim(z) <- c(length(x),length(y))
    if(ssave==TRUE){png(file=file.path("C:","Users","tcvin","OneDrive","Work","DIS","pics","core",paste(nahme,"sd",".png",sep="")))}
    plot.new()
    plot.window(xlim = c(0,1),ylim = c(0,1),asp=NA)
    .filled.contour(x,y, z,levels = seq(0,1.5,0.05), col = terrain.colors(30))
    if(ssave==TRUE){dev.off()}
  }
  if(SK==TRUE){
    xy <- expand.grid(x,y)
    xy2<-cbind(x1value,xy[,1],xy[,2])
    SC<- function(x){(ELKf(x)-f(x))/((varLKf(x))^0.5+0.000000000001)}
    z<-SC(xy2)
    z[which(z> 2)]<-1.99
    dim(z) <- c(length(x),length(y))
    if(ssave==TRUE){png(file=file.path("C:","Users","tcvin","OneDrive","Work","DIS","pics","core",paste(nahme,"sk",".png",sep="")))}
    plot.new()
    plot.window(xlim = c(0,1),ylim = c(0,1),asp=NA)
    .filled.contour(x,y, z,levels = seq(-2.0003,2,0.1), col = colorRampPalette(c('red','chartreuse','navy'),1)(41))
    if(ssave==TRUE){dev.off()}
  }
}

#TD2 takes a slice for fixed x2 value (input x2value) and emulates without any design
TD2<-function(actualplot,expectation,SD,SK,Lev,ssave,nahme,x2value){
  x <- seq(0,1,0.01)
  y <- seq(0,1,0.01)
  xy <- expand.grid(x,y)
  xy2<-cbind(xy[,1],x2value,xy[,2])
  w<-0
  if(actualplot==TRUE){
    z<-f(xy2)
    dim(z) <- c(length(x),length(y))
    
    if(ssave==TRUE){png(file=file.path("C:","Users","tcvin","OneDrive","Work","DIS","pics","core",paste(nahme,"a",".png",sep="")))}
    plot.new()
    plot.window(xlim = c(0,1),ylim = c(0,1),asp=NA)
    .filled.contour(x,y, z,levels = Lev, col = rainbow(25))
    if(ssave==TRUE){dev.off()}
    
  }
  Ef <- B0
  varf <- sigu^2
  r<-function(a){
    exp(-abs(a)^2/theta^2)
  }
  R<-function(a1,a2){
    exp(-abs(a1-a2)^2/theta^2)-exp((-abs(a1)^2-abs(a2)^2)/theta^2)
  }
  nab<-function(x){
    f(x)-Ef
  }
  ELKf <- function(x){
    a<-x[,1]-w
    b<-x[,2]-w
    c<-x[,3]-w
    xK<-cbind(w,x[,2],x[,3])
    xL<-cbind(x[,1],w,w)
    xLK<-cbind(w,xL[,2],w)
    Ef+r(a)*nab(xK)+r(b)*r(c)*nab(xL)-r(a)*r(b)*r(c)*nab(xLK)
  }
  varLKf <- function(x){
    a<-x[,1]-w
    b<-x[,2]-w
    c<-x[,3]-w
    varf*(1-r(a)^2)*(1-r(b)^2*r(c)^2)
  }
  if(expectation==TRUE){
    z<-ELKf(xy2)
    dim(z) <- c(length(x),length(y))
    if(ssave==TRUE){png(file=file.path("C:","Users","tcvin","OneDrive","Work","DIS","pics","core",paste(nahme,"e",".png",sep="")))}
    plot.new()
    plot.window(xlim = c(0,1),ylim = c(0,1),asp=NA)
    .filled.contour(x,y, z,levels = Lev, col = rainbow(25))
    if(ssave==TRUE){dev.off()}
  }
  if(SD==TRUE){
    z<-varLKf(xy2)^0.5
    dim(z) <- c(length(x),length(y))
    if(ssave==TRUE){png(file=file.path("C:","Users","tcvin","OneDrive","Work","DIS","pics","core",paste(nahme,"sd",".png",sep="")))}
    plot.new()
    plot.window(xlim = c(0,1),ylim = c(0,1),asp=NA)
    .filled.contour(x,y, z,levels = seq(0,1.5,0.05), col = terrain.colors(30))
    if(ssave==TRUE){dev.off()}
  }
  if(SK==TRUE){
    xy <- expand.grid(x,y)
    xy2<-cbind(xy[,1],x2value,xy[,2])
    SC<- function(x){(ELKf(x)-f(x))/((varLKf(x))^0.5+0.000000000001)}
    z<-SC(xy2)
    z[which(z> 2)]<-1.99
    dim(z) <- c(length(x),length(y))
    if(ssave==TRUE){png(file=file.path("C:","Users","tcvin","OneDrive","Work","DIS","pics","core",paste(nahme,"sk",".png",sep="")))}
    
    plot.new()
    plot.window(xlim = c(0,1),ylim = c(0,1),asp=NA)
    .filled.contour(x,y, z,levels = seq(-2.0003,2,0.1), col = colorRampPalette(c('red','chartreuse','navy'),1)(41))
    if(ssave==TRUE){dev.off()}
  }
}

#TDD1 takes a slice for fixed x1 value (input x1value) and emulates with design js
TDD1<-function(actualplot,expectation,SD,SK,Lev,ssave,nahme,x1value){
  x <- seq(0,1,0.01)
  y <- seq(0,1,0.01)
  xy <- expand.grid(x,y)
  xy2<-cbind(x1value,xy[,1],xy[,2])
  w<-0
  if(actualplot==TRUE){
    z<-f(xy2)
    dim(z) <- c(length(x),length(y))
    
    if(ssave==TRUE){pdf(file=file.path("C:","Users","tcvin","OneDrive","Work","DIS","pics","core",paste(nahme,"a",".pdf",sep="")))}
    plot.new()
    plot.window(xlim = c(0,1),ylim = c(0,1),asp=NA)
    .filled.contour(x,y, z,levels = Lev, col = rainbow(25))
    js1<-js[which(abs(js[,1]-x1value)<0.1),]
    js2<-js[which(abs(js[,1]-x1value)>0.1 & abs(js[,1]-x1value)<0.2),]
    js3<-js[which(abs(js[,1]-x1value)>0.2 & abs(js[,1]-x1value)<0.3),]
    if(is.matrix(js1)==FALSE){js1<-as.matrix(t(js1))}
    if(is.matrix(js2)==FALSE){js2<-as.matrix(t(js2))}
    if(is.matrix(js3)==FALSE){js3<-as.matrix(t(js3))}
    points(js1[,2],js1[,3],lwd=6)
    points(js2[,2],js2[,3],lwd=3,pch=2)
    points(js3[,2],js3[,3],lwd=3,pch=3)
    if(ssave==TRUE){dev.off()}
    
  }
  x <- seq(0,1,0.01)
  y <- seq(0,1,0.01)
  xy <- expand.grid(x,y)
  xy2<-cbind(x1value,xy[,1],xy[,2])
  Ef <- B0
  varf <- sigu^2
  r<-function(a){
    exp(-abs(a)^2/theta^2)
  }
  R<-function(a1,a2){
    exp(-abs(a1-a2)^2/theta^2)-exp((-abs(a1)^2-abs(a2)^2)/theta^2)
  }
  nab<-function(x){
    f(x)-Ef
  }
  ELKf <- function(x){
    a<-x[,1]-w
    b<-x[,2]-w
    c<-x[,3]-w
    xK<-cbind(w,x[,2],x[,3])
    xL<-cbind(x[,1],w,w)
    xLK<-cbind(w,xL[,2],w)
    Ef+r(a)*nab(xK)+r(b)*r(c)*nab(xL)-r(a)*r(b)*r(c)*nab(xLK)
  }
  varLKf <- function(x){
    a<-x[,1]-w
    b<-x[,2]-w
    c<-x[,3]-w
    varf*(1-r(a)^2)*(1-r(b)^2*r(c)^2)
  }
  covLKfD <- function(x1){
    a1<-x1[,1]-w
    a2<-js[,1]-w
    b1<-x1[,2]-w
    b2<-js[,2]-w
    c1<-x1[,3]-w
    c2<-js[,3]-w
    aa<-expand.grid(a2,a1)
    bb<-expand.grid(b2,b1)
    cc<-expand.grid(c2,c1)
    RRa<-R(aa[,2],aa[,1])
    dim(RRa)<-c(length(js[,2]),length(x1[,2]))
    RRa<-t(RRa)
    rb12<-r(bb[,2])*r(bb[,1])
    dim(rb12)<-c(length(js[,2]),length(x1[,2]))
    rb12<-t(rb12)
    rc12<-r(cc[,2])*r(cc[,1])
    dim(rc12)<-c(length(js[,2]),length(x1[,2]))
    rc12<-t(rc12)
    rbb<-r(bb[,2]-bb[,1])
    dim(rbb)<-c(length(js[,2]),length(x1[,2]))
    rbb<-t(rbb)
    rcc<-r(cc[,2]-cc[,1])
    dim(rcc)<-c(length(js[,2]),length(x1[,2]))
    rcc<-t(rcc)
    RRb<-R(bb[,2],bb[,1])
    dim(RRb)<-c(length(js[,2]),length(x1[,2]))
    RRb<-t(RRb)
    varf*RRa*(rbb*rcc-(rb12*rc12))
  }
  D <- f(js)
  varLKD <- covLKfD(js)
  invvarLKD <- solve(varLKD)
  EDLKf <- function(x){
    A<-t(covLKfD(x))
    C<-(matrix(data=rep(D-ELKf(js),dim(x)[1]),ncol=dim(x)[1]))
    B<-invvarLKD%*%C
    ELKf(x)+colSums(A * B)
  }
  varDLKf <- function(x){
    A<-t(covLKfD(x))
    B<-invvarLKD%*%A
    varLKf(x)-colSums(A * B)
  }
  if(expectation==TRUE){
    z<-EDLKf(xy2)
    dim(z) <- c(length(x),length(y))
    if(ssave==TRUE){pdf(file=file.path("C:","Users","tcvin","OneDrive","Work","DIS","pics","core",paste(nahme,"e",".pdf",sep="")))}
    plot.new()
    plot.window(xlim = c(0,1),ylim = c(0,1),asp=NA)
    .filled.contour(x,y, z,levels = Lev, col = rainbow(25))
    js1<-js[which(abs(js[,1]-x1value)<0.1),]
    js2<-js[which(abs(js[,1]-x1value)>0.1 & abs(js[,1]-x1value)<0.2),]
    js3<-js[which(abs(js[,1]-x1value)>0.2 & abs(js[,1]-x1value)<0.3),]
    if(is.matrix(js1)==FALSE){js1<-as.matrix(t(js1))}
    if(is.matrix(js2)==FALSE){js2<-as.matrix(t(js2))}
    if(is.matrix(js3)==FALSE){js3<-as.matrix(t(js3))}
    points(js1[,2],js1[,3],lwd=6)
    points(js2[,2],js2[,3],lwd=3,pch=2)
    points(js3[,2],js3[,3],lwd=3,pch=3)
    if(ssave==TRUE){dev.off()}
  }
  if(SD==TRUE){
    z<-varDLKf(xy2)^0.5
    dim(z) <- c(length(x),length(y))
    if(ssave==TRUE){pdf(file=file.path("C:","Users","tcvin","OneDrive","Work","DIS","pics","core",paste(nahme,"sd",".pdf",sep="")))}
    plot.new()
    plot.window(xlim = c(0,1),ylim = c(0,1),asp=NA)
    .filled.contour(x,y, z,levels = seq(0,1.5,0.05), col = terrain.colors(30))
    js1<-js[which(abs(js[,1]-x1value)<0.1),]
    js2<-js[which(abs(js[,1]-x1value)>0.1 & abs(js[,1]-x1value)<0.2),]
    js3<-js[which(abs(js[,1]-x1value)>0.2 & abs(js[,1]-x1value)<0.3),]
    if(is.matrix(js1)==FALSE){js1<-as.matrix(t(js1))}
    if(is.matrix(js2)==FALSE){js2<-as.matrix(t(js2))}
    if(is.matrix(js3)==FALSE){js3<-as.matrix(t(js3))}
    points(js1[,2],js1[,3],lwd=6)
    points(js2[,2],js2[,3],lwd=3,pch=2)
    points(js3[,2],js3[,3],lwd=3,pch=3)
    if(ssave==TRUE){dev.off()}
  }
  if(SK==TRUE){
    xy <- expand.grid(x,y)
    xy2<-cbind(x1value,xy[,1],xy[,2])
    SC<- function(x){(EDLKf(x)-f(x))/((varDLKf(x))^0.5+0.000000000001)}
    z<-SC(xy2)
    z[which(z> 2)]<-1.99
    z[which(z< -2)]<-1*-1.99
    dim(z) <- c(length(x),length(y))
    if(ssave==TRUE){pdf(file=file.path("C:","Users","tcvin","OneDrive","Work","DIS","pics","core",paste(nahme,"sk",".pdf",sep="")))}
    
    plot.new()
    plot.window(xlim = c(0,1),ylim = c(0,1),asp=NA)
    .filled.contour(x,y, z,levels = seq(-2.0003,2,0.1), col = colorRampPalette(c('red','chartreuse','navy'),1)(41))
    js1<-js[which(abs(js[,1]-x1value)<0.1),]
    js2<-js[which(abs(js[,1]-x1value)>0.1 & abs(js[,1]-x1value)<0.2),]
    js3<-js[which(abs(js[,1]-x1value)>0.2 & abs(js[,1]-x1value)<0.3),]
    if(is.matrix(js1)==FALSE){js1<-as.matrix(t(js1))}
    if(is.matrix(js2)==FALSE){js2<-as.matrix(t(js2))}
    if(is.matrix(js3)==FALSE){js3<-as.matrix(t(js3))}
    points(js1[,2],js1[,3],lwd=6)
    points(js2[,2],js2[,3],lwd=3,pch=2)
    points(js3[,2],js3[,3],lwd=3,pch=3)
    if(ssave==TRUE){dev.off()}
  }
}

#TDD2 takes a slice for fixed x2 value (input x1value) and emulates with design js
TDD2<-function(actualplot,expectation,SD,SK,Lev,ssave,nahme,x2value){
  x <- seq(0,1,0.01)
  y <- seq(0,1,0.01)
  xy <- expand.grid(x,y)
  xy2<-cbind(xy[,1],x2value,xy[,2])
  w<-0
  if(actualplot==TRUE){
    z<-f(xy2)
    dim(z) <- c(length(x),length(y))
    
    if(ssave==TRUE){pdf(file=file.path("C:","Users","tcvin","OneDrive","Work","DIS","pics","core",paste(nahme,"a",".pdf",sep="")))}
    plot.new()
    plot.window(xlim = c(0,1),ylim = c(0,1),asp=NA)
    .filled.contour(x,y, z,levels = Lev, col = rainbow(25))
    js1<-js[which(abs(js[,2]-x2value)<0.1),]
    js2<-js[which(abs(js[,2]-x2value)>0.1 & abs(js[,2]-x2value)<0.2),]
    js3<-js[which(abs(js[,2]-x2value)>0.2 & abs(js[,2]-x2value)<0.3),]
    if(is.matrix(js1)==FALSE){js1<-as.matrix(t(js1))}
    if(is.matrix(js2)==FALSE){js2<-as.matrix(t(js2))}
    if(is.matrix(js3)==FALSE){js3<-as.matrix(t(js3))}
    points(js1[,1],js1[,3],lwd=6)
    points(js2[,1],js2[,3],lwd=3,pch=2)
    points(js3[,1],js3[,3],lwd=3,pch=3)
    if(ssave==TRUE){dev.off()}
    
  }
  x <- seq(0,1,0.01)
  y <- seq(0,1,0.01)
  xy <- expand.grid(x,y)
  xy2<-cbind(xy[,1],x2value,xy[,2])
  Ef <- B0
  varf <- sigu^2
  r<-function(a){
    exp(-abs(a)^2/theta^2)
  }
  R<-function(a1,a2){
    exp(-abs(a1-a2)^2/theta^2)-exp((-abs(a1)^2-abs(a2)^2)/theta^2)
  }
  nab<-function(x){
    f(x)-Ef
  }
  ELKf <- function(x){
    a<-x[,1]-w
    b<-x[,2]-w
    c<-x[,3]-w
    xK<-cbind(w,x[,2],x[,3])
    xL<-cbind(x[,1],w,w)
    xLK<-cbind(w,xL[,2],w)
    Ef+r(a)*nab(xK)+r(b)*r(c)*nab(xL)-r(a)*r(b)*r(c)*nab(xLK)
  }
  varLKf <- function(x){
    a<-x[,1]-w
    b<-x[,2]-w
    c<-x[,3]-w
    varf*(1-r(a)^2)*(1-r(b)^2*r(c)^2)
  }
  covLKfD <- function(x1){
    a1<-x1[,1]-w
    a2<-js[,1]-w
    b1<-x1[,2]-w
    b2<-js[,2]-w
    c1<-x1[,3]-w
    c2<-js[,3]-w
    aa<-expand.grid(a2,a1)
    bb<-expand.grid(b2,b1)
    cc<-expand.grid(c2,c1)
    RRa<-R(aa[,2],aa[,1])
    dim(RRa)<-c(length(js[,2]),length(x1[,2]))
    RRa<-t(RRa)
    rb12<-r(bb[,2])*r(bb[,1])
    dim(rb12)<-c(length(js[,2]),length(x1[,2]))
    rb12<-t(rb12)
    rc12<-r(cc[,2])*r(cc[,1])
    dim(rc12)<-c(length(js[,2]),length(x1[,2]))
    rc12<-t(rc12)
    rbb<-r(bb[,2]-bb[,1])
    dim(rbb)<-c(length(js[,2]),length(x1[,2]))
    rbb<-t(rbb)
    rcc<-r(cc[,2]-cc[,1])
    dim(rcc)<-c(length(js[,2]),length(x1[,2]))
    rcc<-t(rcc)
    RRb<-R(bb[,2],bb[,1])
    dim(RRb)<-c(length(js[,2]),length(x1[,2]))
    RRb<-t(RRb)
    varf*RRa*(rbb*rcc-(rb12*rc12))
  }
  D <- f(js)
  varLKD <- covLKfD(js)
  invvarLKD <- solve(varLKD)
  EDLKf <- function(x){
    A<-t(covLKfD(x))
    C<-(matrix(data=rep(D-ELKf(js),dim(x)[1]),ncol=dim(x)[1]))
    B<-invvarLKD%*%C
    ELKf(x)+colSums(A * B)
  }
  varDLKf <- function(x){
    A<-t(covLKfD(x))
    B<-invvarLKD%*%A
    varLKf(x)-colSums(A * B)
  }
  if(expectation==TRUE){
    z<-EDLKf(xy2)
    dim(z) <- c(length(x),length(y))
    if(ssave==TRUE){pdf(file=file.path("C:","Users","tcvin","OneDrive","Work","DIS","pics","core",paste(nahme,"e",".pdf",sep="")))}
    plot.new()
    plot.window(xlim = c(0,1),ylim = c(0,1),asp=NA)
    .filled.contour(x,y, z,levels = Lev, col = rainbow(25))
    js1<-js[which(abs(js[,2]-x2value)<0.1),]
    js2<-js[which(abs(js[,2]-x2value)>0.1 & abs(js[,2]-x2value)<0.2),]
    js3<-js[which(abs(js[,2]-x2value)>0.2 & abs(js[,2]-x2value)<0.3),]
    if(is.matrix(js1)==FALSE){js1<-as.matrix(t(js1))}
    if(is.matrix(js2)==FALSE){js2<-as.matrix(t(js2))}
    if(is.matrix(js3)==FALSE){js3<-as.matrix(t(js3))}
    points(js1[,1],js1[,3],lwd=6)
    points(js2[,1],js2[,3],lwd=3,pch=2)
    points(js3[,1],js3[,3],lwd=3,pch=3)
    if(ssave==TRUE){dev.off()}
  }
  if(SD==TRUE){
    z<-varDLKf(xy2)^0.5
    dim(z) <- c(length(x),length(y))
    if(ssave==TRUE){pdf(file=file.path("C:","Users","tcvin","OneDrive","Work","DIS","pics","core",paste(nahme,"sd",".pdf",sep="")))}
    plot.new()
    plot.window(xlim = c(0,1),ylim = c(0,1),asp=NA)
    .filled.contour(x,y, z,levels = seq(0,1.5,0.05), col = terrain.colors(30))
    js1<-js[which(abs(js[,2]-x2value)<0.1),]
    js2<-js[which(abs(js[,2]-x2value)>0.1 & abs(js[,2]-x2value)<0.2),]
    js3<-js[which(abs(js[,2]-x2value)>0.2 & abs(js[,2]-x2value)<0.3),]
    if(is.matrix(js1)==FALSE){js1<-as.matrix(t(js1))}
    if(is.matrix(js2)==FALSE){js2<-as.matrix(t(js2))}
    if(is.matrix(js3)==FALSE){js3<-as.matrix(t(js3))}
    points(js1[,1],js1[,3],lwd=6)
    points(js2[,1],js2[,3],lwd=3,pch=2)
    points(js3[,1],js3[,3],lwd=3,pch=3)
    if(ssave==TRUE){dev.off()}
  }
  if(SK==TRUE){
    xy <- expand.grid(x,y)
    xy2<-cbind(xy[,1],x2value,xy[,2])
    SC<- function(x){(EDLKf(x)-f(x))/((varDLKf(x))^0.5+0.000000000001)}
    z<-SC(xy2)
    z[which(z> 2)]<-1.99
    z[which(z< -2)]<-1*-1.99
    dim(z) <- c(length(x),length(y))
    if(ssave==TRUE){pdf(file=file.path("C:","Users","tcvin","OneDrive","Work","DIS","pics","core",paste(nahme,"sk",".pdf",sep="")))}
    
    plot.new()
    plot.window(xlim = c(0,1),ylim = c(0,1),asp=NA)
    .filled.contour(x,y, z,levels = seq(-2.0003,2,0.1), col = colorRampPalette(c('red','chartreuse','navy'),1)(41))
    js1<-js[which(abs(js[,2]-x2value)<0.1),]
    js2<-js[which(abs(js[,2]-x2value)>0.1 & abs(js[,2]-x2value)<0.2),]
    js3<-js[which(abs(js[,2]-x2value)>0.2 & abs(js[,2]-x2value)<0.3),]
    if(is.matrix(js1)==FALSE){js1<-as.matrix(t(js1))}
    if(is.matrix(js2)==FALSE){js2<-as.matrix(t(js2))}
    if(is.matrix(js3)==FALSE){js3<-as.matrix(t(js3))}
    points(js1[,1],js1[,3],lwd=6)
    points(js2[,1],js2[,3],lwd=3,pch=2)
    points(js3[,1],js3[,3],lwd=3,pch=3)
    if(ssave==TRUE){dev.off()}
  }
}

##########
# FIGURES
##########
# As the report required the true Aradopsis Model we cannot put that here to see the same plots. However we can example
# the code with the following toy function:
f <- function(x){tan(x[,1]*x[,3])+sin(x[,2]*x[,3])*log(x[,1]+2)+cos(3.8*x[,1])+0.01}
B0 = 0.34
sigu = 1.5
theta = 0.8
Lev=seq(-1,1.1,0.1)
library(lhs)
js<-maximinLHS(16,3) #The LHS used in the report was maximined by easy manipulation of code in figure 3.4 for 3D

# The values x1value and x2value can be change to any number from 0 to 1 to view emulation through the entire plane
TD1(actualplot=TRUE,expectation=TRUE,SD=TRUE,SK=TRUE,Lev=seq(-1,1.7,0.1),ssave=ss,nahme='three',x1value=0)
TD2(actualplot=TRUE,expectation=TRUE,SD=TRUE,SK=FALSE,Lev=seq(-1,1.7,0.1),ssave=ss,nahme='three',x2value=0)

TDD1(actualplot=TRUE,expectation=TRUE,SD=TRUE,SK=FALSE,Lev=seq(-1,1.7,0.1),ssave=ss,nahme='three',x1value=0)
TDD2(actualplot=TRUE,expectation=TRUE,SD=TRUE,SK=FALSE,Lev=seq(-1,1.7,0.1),ssave=ss,nahme='three',x2value=0)

