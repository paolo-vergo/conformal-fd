## fData #############################?

N = 20
P = 1e2
grid = seq( 0, 1, length.out = P )
C = roahd::exp_cov_function( grid, alpha = 0.3, beta = 0.4 )
values = roahd::generate_gauss_fdata( N,
                               centerline = sin( 2 * pi * grid ),
                               Cov = C )
fD = roahd::fData( grid, values )
x0=list(as.list(grid))
fun=mean_lists()
final.fData = conformal.fun.split(NULL,NULL, fD, x0, fun$train.fun, fun$predict.fun,
                             alpha=0.1,
                             split=NULL, seed=FALSE, randomized=FALSE,seed_tau=FALSE,
                             verbose=TRUE, training_size=0.5,s_type="alpha-max")
plot_fun(final.fData)

###  mfData ###################################

N = 1e2
P = 1e3
t0 = 0
t1 = 1
grid = seq( t0, t1, length.out = P )
C = roahd::exp_cov_function( grid, alpha = 0.3, beta = 0.4 )
Data_1 = roahd::generate_gauss_fdata( N, centerline = sin( 2 * pi * grid ), Cov = C )
Data_2 = roahd::generate_gauss_fdata( N, centerline = log(1+ 2 * pi * grid ), Cov = C )
mfD=roahd::mfData( grid, list( Data_1, Data_2 ) )
x0=list(as.list(grid))
fun=mean_lists()
final.mfData = conformal.fun.split(NULL,NULL, mfD, x0, fun$train.fun, fun$predict.fun,
                             alpha=0.1,
                             split=NULL, seed=FALSE, randomized=FALSE,seed_tau=FALSE,
                             verbose=TRUE, training_size=0.5,s_type="alpha-max")
h=plot_fun(final.mfData)

### fd ###########################################

daybasis <- fda::create.fourier.basis(c(0, 365), nbasis=65)
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
final.fd = conformal.fun.split(NULL,fda::day.5, tempfd, x0, fun$train.fun, fun$predict.fun,
                             alpha=0.1,
                             split=NULL, seed=FALSE, randomized=FALSE,seed_tau=FALSE,
                             verbose=TRUE, training_size=0.5,s_type="alpha-max")
plot_fun(final.fd)

