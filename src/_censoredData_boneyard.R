## censored data
tmp=subset(PAH.dat.melt,variable=="naphthalene")
library(NADA)
#  (Kaplan-Meier method)
with(tmp,Cen(value,cen))
test.NADA=cenfit(tmp$value, tmp$cen)

## non NADA function
Surv=Surv(tmp$value, !tmp$cen, type="left")
#  (ROS)
myros=ros(tmp$value, tmp$cen)
plot(myros)
summary(myros)
mean(myros); sd(myros)
quantile(myros); median(myros)
as.data.frame(myros)

# from ROS code
pp = hc.ppoints(tmp$value, tmp$cen, "na.omit");# arbitrarly assigns pp values to censored data
pp.nq = qnorm(pp[!tmp$cen])
obs.transformed = log(tmp$value[!tmp$cen])
hc = lm(obs.transformed ~ pp.nq, na.action="na.omit")
oldClass(hc) = c("ros", "lm")
hc$obs      = tmp$value
hc$modeled  = tmp$value 
hc$pp       = pp
hc$censored = tmp$cen 
hc$reverseT = "exp"
hc$forwardT = "log"
predict(hc, data.frame(pp.nq=qnorm(pp[tmp$cen])), na.action="na.omit")


obs=tmp$value
censored=tmp$cen
# cohn function
uncen = obs[!censored]
cen = obs[censored]
A = B = C = P = numeric()
limit = sort(unique(cen))
a = length(uncen[uncen < limit[1]])
if (a > 0) {
  limit = c(0, limit)
}
i = length(limit)
A[i] = length(uncen[uncen >= limit[i]])
B[i] = length(obs[obs <= limit[i]]) - length(uncen[uncen == 
                                                     limit[i]])
C[i] = length(cen[cen == limit[i]])
P[i] = A[i]/(A[i] + B[i])
i = i - 1
while (i > 0) {
  A[i] = length(uncen[uncen >= limit[i] & uncen < limit[i + 
                                                          1]])
  B[i] = length(obs[obs <= limit[i]]) - length(uncen[uncen == 
                                                       limit[i]])
  C[i] = length(cen[cen == limit[i]])
  P[i] = P[i + 1] + ((A[i]/(A[i] + B[i])) * (1 - P[i + 
                                                     1]))
  i = i - 1
}
return(list(A = A, B = B, C = C, P = P, limit = limit))


hc.ppoints.uncen
hc.ppoints.cen

cn=cohn(obs,censored)
hc.ppoints.uncen(obs,censored,cn,"na.omit")
hc.ppoints.cen(obs,censored,cn,"na.omit")

pp = numeric(length(obs))
pp[!censored]=hc.ppoints.uncen(obs,censored,cn,"na.omit")
