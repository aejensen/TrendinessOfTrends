#Data obtained from "fra Nøgletal - Danskernes rygevaner 2018"
#at https://www.sst.dk/da/udgivelser/2019/danskernes-rygevaner-2018
year <- c(1998, 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018)
p <- c(34.6, 34.1, 33.5, 32.3, 31.0, 30.0, 27.1, 28.0, 27.7, 28.5, 28.0, 24.3, 23.4, 22.3, 22.6, 21.0, 22.5, 21.1, 21.6, 23.1)
smoking <- data.frame(year = year, p = p)

plot(p ~ year, smoking, pch = 19)
save(smoking, file="smoking.RData")
