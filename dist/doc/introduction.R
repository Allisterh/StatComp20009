## -----------------------------------------------------------------------------
library(StatComp20009)
data(dt.d)
rt.log<-log(dt.d[2:dim(dt.d)[1],]/dt.d[1:(dim(dt.d)[1]-1),])
r.B1.day<-getnfac(rt.log,20,'PC1')$ic
r.B2.day<-getnfac(rt.log,20,'PC2')$ic
r.B3.day<-getnfac(rt.log,20,'PC3')$ic
r.Y.day<-nof(rt.log,3)
knitr::kable(cbind(r.B1.day,r.B2.day,r.B3.day,r.Y.day))

## -----------------------------------------------------------------------------
L<-30
k<-8
r.B1<-r.B2<-r.B3<-r.Y<-numeric(k)
for (i in 1:k) {
  L<-30
  rt.log<-log(dt.d[(2:L)+L*(i-1),]/dt.d[(1:(L-1))+L*(i-1),])
  L<-L-1
  r.B1[i]<-getnfac(rt.log,20,'PC1')$ic
  r.B2[i]<-getnfac(rt.log,20,'PC2')$ic
  r.B3[i]<-getnfac(rt.log,20,'PC3')$ic
  r.Y[i]<-nof(rt.log,3)
}
knitr::kable(cbind(r.B1,r.B2,r.B3,r.Y))

