# Math

2+3
10-5
3-7
2*5
5/8
10/2
2^2

x <- 5
y <- 10
x/y
y*x
x+2

# character

"ali"
"emam"
ali
"ali" + "emam"
"ali" - "emam"
"5" + "5"
5 + "5"
x <- "ali"
y <- "emam"

# vector

c(1,2,3,4,5)
c("ali", "emam")
c(1,2,3, "ali") # only one type of data

x <- c(1,2,3,4,5)
y <- c("ali", "emam")
xx <- c(10,20,30,40,50)
xx + x

# functions

sqrt(25)
sqrt("25")

mean(x)
mean(y)
mean(xx)
xx <- c(10,20,30,40,50, NA)
mean(xx)
mean(xx, na.rm = T)
median(xx, na.rm = T)
summary(xx)
# def. new function
xx <- c(10,20,30,40,40, 40,50)
yy <- table(xx)
yy
which.max(yy)

get.mode <- function(x){
  y <- table(x)
  mode <- which.max(y)
  return(mode)
}

get.mode(xx)

# data.frame

gene1 <- c(20,30,40,60)
gene2 <- c(200,350,40,30)
gene3 <- c(56,156,23,75)
gene4 <- c(30,11,19,56)
df.genes <- data.frame(g1 = gene1, g2= gene2, g3= gene3, g4 = gene4)
View(df.genes)
mean(df.genes$g2)
mean(df.genes$g3)

df.genes$g1 + df.genes$g4

apply(df.genes, 2, mean) # col. by genes

apply(df.genes, 1, mean) # row. by sample

plot(df.genes)
hist(df.genes$g1)
hist(df.genes$g1,freq = F)
#install.packages("ggplot2")

library(ggplot2)
ggplot(df.genes,aes(x = g2))+ geom_histogram(bins = 5)
