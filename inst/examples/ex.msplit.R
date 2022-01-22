## fData #############################

N = 30
P = 60
grid = seq( 0, 1, length.out = P )
C = roahd::exp_cov_function( grid, alpha = 0.3, beta = 0.4 )
values = roahd::generate_gauss_fdata( N,
                                      centerline = sin( 2 * pi * grid ),
                                      Cov = C )
fD = roahd::fData( grid, values )
x0=list(as.list(grid))
fun=mean_lists()
B=4

final.multi=conformal.fun.msplit(NULL,NULL, fD, x0, fun$train.fun, fun$predict.fun,
                                 alpha=0.1,
                                 split=NULL, seed=FALSE, randomized=FALSE,seed_beta=FALSE,
                                 verbose=FALSE, training_size=NULL,s_type="st-dev",B,lambda=0,
                                 tau = 0)

###  mfData ###################################

N = 50
P = 50
t0 = 0
t1 = 1
grid = seq( t0, t1, length.out = P )
C = roahd::exp_cov_function( grid, alpha = 0.3, beta = 0.4 )
Data_1 = roahd::generate_gauss_fdata( N, centerline = sin( 2 * pi * grid ), Cov = C )
Data_2 = roahd::generate_gauss_fdata( N, centerline = log(1+ 2 * pi * grid ), Cov = C )
mfD=roahd::mfData( grid, list( Data_1, Data_2 ) )
x0=list(as.list(grid))
fun=mean_lists()
B=4
final.multi=conformal.fun.msplit(NULL,NULL, mfD, x0, fun$train.fun, fun$predict.fun,
                                 alpha=0.1,
                                 split=NULL, seed=FALSE, randomized=FALSE,seed_beta=FALSE,
                                 verbose=FALSE, training_size=NULL,s_type="st-dev",B,lambda=0,
                                 tau = 1-(B+1)/(2*B))


### fd ###########################################
daybasis <- fda::create.fourier.basis(c(0, 365), nbasis=25)
tempfd <- fda::smooth.basis(fda::day.5, fda::CanadianWeather$dailyAv[,,"Temperature.C"],daybasis)$fd
Lbasis <- fda::create.constant.basis(c(0, 365))
Lcoef <- matrix(c(0,(2*pi/365)^2,0),1,3)
bfdobj <- fda::fd(Lcoef,Lbasis)
bwtlist <- fda::fd2list(bfdobj)
harmaccelLfd <- fda::Lfd(3, bwtlist)
Ltempmat <- fda::eval.fd(fda::day.5, tempfd, harmaccelLfd)
t=1:365
x0=list(as.list(grid))
fun=mean_lists()
B=4
final.multi=conformal.fun.msplit(NULL,fda::day.5, tempfd, x0, fun$train.fun, fun$predict.fun,
                        alpha=0.1,
                        split=NULL, seed=FALSE, randomized=FALSE,seed_beta=FALSE,
                        verbose=FALSE, training_size=NULL,s_type="st-dev",B,lambda=0,
                        tau = 0.2)


