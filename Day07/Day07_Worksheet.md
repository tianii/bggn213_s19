Day07:Functions and packages
================
Tiani Louis
4/24/2019

More on function writing
------------------------


``` r
source("http://tinyurl.com/rescale-R")
```

Test the *rescale()* function

``` r
rescale(1:10)
```

    ##  [1] 0.0000000 0.1111111 0.2222222 0.3333333 0.4444444 0.5555556 0.6666667
    ##  [8] 0.7777778 0.8888889 1.0000000

``` r
x <- c(1:10, "string")
!is.numeric(x)
```

    ## [1] TRUE

Function practice
-----------------

Write a function to identify NA elements in two vectors Start with a simple example where I know what the answer should be

``` r
x <- c( 1, 2, NA, 3, NA)
y <- c(NA, 3, NA, 3, 4)
```

``` r
is.na(x)
```

    ## [1] FALSE FALSE  TRUE FALSE  TRUE

``` r
is.na(y)
```

    ## [1]  TRUE FALSE  TRUE FALSE FALSE

this will tell you when both x and y are NA
===========================================

``` r
is.na(x) & is.na(y)
```

    ## [1] FALSE FALSE  TRUE FALSE FALSE

Take the sum to find how many

``` r
sum( is.na(x) & is.na(y))
```

    ## [1] 1

Take this and use it as the body of the function

``` r
both_na <- function(x,y) {
if(length(x) != length(y)) {
stop("UHHHH, They really should be the same length")
 }
sum( is.na(x) & is.na(y))  
}
```

TEST IT!

``` r
both_na(x,y)
```

    ## [1] 1

``` r
both_na(c(NA,NA,NA), c(1,NA,NA))
```

    ## [1] 2

These are different lengths and this code should throw an error &gt; R will recycle the length of the first vector to produce an answer

``` r
# both_na(c(NA,NA,NA), c(1,NA,NA,NA))
```

Check the length of the vectors using length()

``` r
length(x) == length(y)
```

    ## [1] TRUE

Add some cute extras!

``` r
x <- c( 1, 2, NA, 3, NA)
y <- c(NA, 3, NA, 3, 4)

both_na <- function(x,y) {
  if(length(x) != length(y)) {
 stop("UHHHH, They really should be the same length")
 }
sum( is.na(x) & is.na(y))  
}

both_na(x,y)
```

    ## [1] 1

``` r
# student 1
s_1 <- c(100,100,100,100,100,100,100,90)

#student 2
s_2 <- c(100, NA, 90, 90,90,90, 97, 80)
```

``` r
grade <- function(x) {
  min_g <- min(x)
 sum_g <- sum(x) - min_g
 average <- sum_g/(length(x)-1)
return(average)
}
```

``` r
grade(s_1)
```

    ## [1] 100

``` r
grade(s_2)
```

    ## [1] NA

``` r
grade2 <- function(x) {
  min_g <- min(x)
 average <- mean(x, na.rm =T)
return(average)
}
```

``` r
grade2(s_2)
```

    ## [1] 91

``` r
grade3 <- function(x) {
  min_g <- min(x)
  sum_g <- sum(x) - min_g
  average <- sum_g/(length(x)-1)
  if (is.na(x)) {
  min_g <- min(x)
 average2 <- mean(x, na.rm =T)
return(average2)
  }
 
}
```

``` r
grade3(s_2)
```

    ## Warning in if (is.na(x)) {: the condition has length > 1 and only the first
    ## element will be used

``` r
grade <-function(x) {
  (sum(x, na.rm=T) - min(x, na.rm=T))/(length(x)-1)
}
```

``` r
grade(s_2)
```

    ## [1] 79.57143

``` r
# student_homework <- read.csv("~/Desktop/BGGN213/Day07/student_homework.csv")
```

``` r
url <- "https://tinyurl.com/gradeinput"
students <- read.csv(url, row.names = 1)
head(students)
```

    ##           hw1 hw2 hw3 hw4 hw5
    ## student-1 100  73 100  88  79
    ## student-2  85  64  78  89  78
    ## student-3  83  69  77 100  77
    ## student-4  88  NA  73 100  76
    ## student-5  88 100  75  86  79
    ## student-6  89  78 100  89  77

``` r
grade(students[5,])
```

    ## [1] 88.25

``` r
ans <- apply(students, 1, grade)
```

``` r
sort(ans, decreasing =T)
```

    ##  student-7  student-8 student-13  student-1 student-12 student-16 
    ##      94.00      93.75      92.25      91.75      91.75      89.50 
    ##  student-6  student-5 student-17  student-9 student-14 student-11 
    ##      89.00      88.25      88.00      87.75      87.75      86.00 
    ##  student-3 student-19 student-20  student-2 student-18  student-4 
    ##      84.25      82.75      82.75      82.50      72.75      66.00 
    ## student-15 student-10 
    ##      62.50      61.00

The function to end all functions find the intersection of two samples

``` r
x <- df1$IDs
y <- df2$IDs
z <- df3$IDs
intersect(x,y)
```

    ## [1] "gene2" "gene3"

``` r
x
```

    ## [1] "gene1" "gene2" "gene3"

``` r
y
```

    ## [1] "gene2" "gene4" "gene3" "gene5"

``` r
z
```

    ## [1] "gene2" "gene2" "gene5" "gene5"

``` r
 x[x %in% z]
```

    ## [1] "gene2"

column bind cbind() row bind rbind() %in% finds the common elements between two vectors

``` r
gene_intersect <- function(x, y) {
  cbind(x[x %in% y], y[y %in% x])
}
```

``` r
merge(df1, df2, by="IDs")
```

    ##     IDs exp.x exp.y
    ## 1 gene2     1    -2
    ## 2 gene3     1     1

``` r
# ggplot(df1, aes( x = IDs)) + geom_path(aes(y = exp), color="red")
```
