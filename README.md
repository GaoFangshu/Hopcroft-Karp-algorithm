# Hopcroft-Karp-algorithm
Implementation of Hopcroft-Karp algorithm in R

Here is an example:
```
source('C:/Users/Fangshu Gao/Desktop/Advisor/Hopcroft_Karp.R', encoding = 'UTF-8')
test1 <- matrix(c(1,1,0,0,1,0,1,
                  1,0,0,0,0,1,0,
                  0,1,0,0,0,0,0,
                  0,1,0,1,0,0,0,
                  0,1,0,1,0,0,0,
                  0,0,1,1,1,0,0,
                  0,0,0,0,0,1,0), nrow = 7, ncol = 7, byrow = TRUE)
test2 <- matrix(c(1,0,1,
                  0,1,0,
                  1,0,0), nrow = 3, ncol = 3, byrow = TRUE)
MMHK(test1)
> [1] 3 2 4 4 6 3 7 6 2 1 1 5
MMHK(test2)
> [1] 2 2 3 1 1 3
```

The result shows that matrix `test1` has a maximum cardinality matching {(3,2),(4,4),(6,3),(7,6),(2,1),(1,5)} and matrix `test2` has a perfect matching {(2,2),(3,1),(1,3)}.

More details can be found at my [blog](http://gaofangshu.com/blog/2017/04/11/Implementation-of-Hopcroft-Karp-Algorithm/).
