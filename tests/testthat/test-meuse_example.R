data("meuse.all", package = "gstat")

test_that("fitting model with coords specified work", {
  m <- DSSP(
    formula=log(zinc)~1, data=meuse.all,N=1000, function(x) -2*log(1+x),
    fitted.values=TRUE,pars=c(0.001,0.001),coords=~x+y)
  expect_true(class(m)=="list")
})

sp::coordinates(meuse.all) = ~x+y
X<-scale(sp::coordinates(meuse.all))
X.train<-X[1:155,]
Y<-scale(log(meuse.all$zinc))
Y.train<-Y[1:155]
X.pred<-X[156:164,]
N<-10000

center.y<-attr(Y,"scaled:center")
scale.y<-attr(Y,"scaled:scale")

test_that("fitting model work", {
  meuse.fit<<-DSSP(formula=log(zinc)~1, data=meuse.all[1:155,],N=N, function(x) -2*log(1+x),fitted.values=TRUE,pars=c(0.001,0.001))
  ETA<<-meuse.fit$eta
  expect_true(is.numeric(ETA))
  DELTA<<-meuse.fit$delta
  expect_true(is.numeric(DELTA))
  NU<<-meuse.fit$nu
  expect_true(is.numeric(NU))
})

test_that("fitted values are good", {
  Yhat<-rowMeans(exp(NU*scale.y+center.y))
  Ytrue<-meuse.all$zinc[1:155]
  mse<-mean((Ytrue - Yhat)^2)
  expect_true(mse < 40000)
  expect_true(cor(Ytrue, Yhat) > 0.8)
})

test_that("predictions are good", {
  Y.pred<-DSSP.predict(meuse.fit,X.pred)
  Y.pred<-exp(Y.pred*scale.y+center.y)
  Y.pred<-rowMeans(Y.pred)
  Y.true<-meuse.all$zinc[156:164]
  expect_true(mean((Y.true - Y.pred)^2) < 20000)
  expect_true(cor(Y.pred, Y.true) > 0.8)
})
