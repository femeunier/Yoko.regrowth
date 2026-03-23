rm(list = ls())

load("./outputs/Yoko_default.RData")

matplot(datum$year,
        apply(datum$szpft$agb[,2:11,c(2,3,4,18)],c(1,3),sum),
        type = "l",lty = 1,col = c("red","green","blue","black"))

matplot(datum$year,
        datum$szpft$lai[,12,c(2,3,4,18)],
        type = "l",lty = 1,col = c("red","green","blue","black"))

matplot(datum$year,
        datum$szpft$gpp[,12,c(2,3,4,18)],type = "l",lty = 1,col = c("red","green","blue","black"))

plot(datum$year,
     datum$emean$nep,type = 'l')
abline(h = 0,lty = 2,col = "red")


matplot((datum$emean$soil.water[,c(1,16)]),type = "l")
