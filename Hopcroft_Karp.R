# 2017.4.7
# Hopcroft–Karp algorithm


MMHK <- function(m_input) { # maximum matching with Hopcroft-Karp algorithm
  
  M <- m_input
  
  
  ## prepare functions ---------
  
  
  # input: c(a,b) 
  # output: c(a,c),c(a,d)... as matrix
  shift_on_row <- function(x) {
    index_ones <- c(1:length(M[x[1], ]))[M[x[1], ] == 1]
    index_ones <- index_ones[!(index_ones %in% x[2])]
    output <- sapply(index_ones, function(i) c(x[1], i))
    return(output)
  }
  
  # input: c(a,b) 
  # output: c(c,b) in MM
  shift_on_col <- function(x, MM) {
    index_col <- x[2]
    index_col_output <- seq(2, length(MM), 2)[MM[seq(2, length(MM), 2)] %in% index_col]
    if (length(index_col_output) == 0) {
      return(NA)
    } else {
      output <- MM[c(index_col_output - 1, index_col_output)]
      return(output)
    }
  }
  
  # initiate the match set of P
  init_match <- function(){
    P <- c()
    M0 <- M
    rownames(M0) <- c(1:nrow(M))
    colnames(M0) <- c(1:ncol(M))
    for (i in 1:nrow(M)) {
      M0 <- as.matrix(M0)
      if (sum(M0[1, ]) > 0) {
        P <- c(P,c(i, (as.integer(colnames(M0))[M0[1, ] == 1])[1]))
        M0 <- M0[-c(1:nrow(M0))[as.integer(rownames(M0)) %in% P[length(P) - 1]], 
                 -c(1:ncol(M0))[as.integer(colnames(M0)) %in% P[length(P)]]]
      } else {
        if (nrow(M0) == 0) {
          break
        } else {
          zerorow <- c(1:nrow(M0))[as.vector(rowSums(M0)) == 0]
          zerocol <- c(1:ncol(M0))[as.vector(colSums(M0)) == 0]
          M0 <- M0[-zerorow[1], -zerocol[1]]
        }
      }
    }
    return(P)
  }
  
  # set entry of the tree
  # input: MM
  # output: c(a,0),c(b,0)... as matrix
  set_entry <- function(MM) {
    index_row <- c(1:nrow(M))[!c(1:nrow(M)) %in% MM[seq(1, length(MM), 2)]]
    output <- sapply(index_row, function(i) c(i, 0))
    return(output)
  }
  
  # look for path
  # input: start_point c(a,0),c(b,0)... as matrix
  # output: path list
  find_path <- function(x, path_list, MM) { # x is initially start_point
    all_layer_1 <- c()
    all_layer_2 <- c()
    
    # if (length(shift_on_row(x[, i])) > 0) {
    #   # can shift row this time
    #   if (!is.na(shift_on_col(layer_1[, j]))) {
    #     # can shift column this time
    #   } else {
    #     # Augmenting Path is found!
    #   }
    # } else {
    #   # path finished but not a Augmenting Path
    # }
    
    for (i in 1:ncol(x)) {
      ## add start_point[, i] to list
      if (length(shift_on_row(x[, i])) > 0) { # can shift row this time
        layer_1 <- shift_on_row(x[, i])
        all_layer_1 <- c(all_layer_1, as.vector(layer_1))
        for (j in 1:ncol(layer_1)) {
          if (!is.na(shift_on_col(layer_1[, j], MM)[1])) { # can shift column this time
            ## add layer_1[, j] to list
            layer_2 <- shift_on_col(layer_1[, j], MM)
            all_layer_2 <- c(all_layer_2, layer_2) # collect all possible layer_2
          } else { # Augmenting Path is found!
            # pass
          } 
        }
      } else { # path finished but not a Augmenting Path
        # pass
      }
    }
    if (length(all_layer_1) == 0 & length(all_layer_2) == 0) { # all path finished!
      return(path_list)
    } else { # start next iteration
      if (length(all_layer_1) != 0) {
        path_list <- c(path_list, list(all_layer_1))
      }
      if (length(all_layer_2) != 0) {
        path_list <- c(path_list, list(all_layer_2))
        new_start_points <- (matrix(all_layer_2, nrow = 2)) # update start points
        find_path(new_start_points, path_list, MM)
      } else {
        return(path_list)
      }
    }
  }
  
  # get Augmenting Path
  # input: path list
  # output: vector of Augmenting Path
  get_path <- function(x) {
    path <- c()
    find_by_row <- function(x, layer) {
      index_1 <- seq(1, length(layer), 2)[layer[seq(1, length(layer), 2)] %in% x[1]][1]
      output <- layer[c(index_1, index_1 + 1)]
      return(output)
    }
    find_by_col <- function(x, layer) {
      index_1 <- seq(2, length(layer), 2)[layer[seq(2, length(layer), 2)] %in% x[2]][1]
      output <- layer[c(index_1 - 1, index_1)]
      return(output)
    }
    delete_edge_by_col <- function(x, layer) {
      index_1 <- seq(2, length(layer), 2)[x[2] %in% layer[seq(2, length(layer), 2)]][1]
      output <- layer[c(-(index_1 - 1), -index_1)]
      return(output)
    }
    if (length(x) %% 2 == 1) { # list not ended with Augmenting Path or no Augmenting Path
      for (i in 1:(length(x[[length(x)]])/2)) {
        tested_edge <- c(x[[length(x)]][i * 2 - 1], x[[length(x)]][i * 2])
        finded_edge <- find_by_col(tested_edge, x[[length(x) - 1]])
        if (!is.na(finded_edge)[1]) {
          x[[length(x)]] <- x[[length(x)]][c(-(i * 2 - 1), -(i * 2))]
        }
        while (!is.na(finded_edge)[1]) {
          x[[length(x) - 1]] <- delete_edge_by_col(tested_edge, x[[length(x) - 1]])
          if(length(x[[length(x) - 1]]) == 0) {
            finded_edge <- NA
          } else {
            finded_edge <- find_by_col(tested_edge, x[[length(x) - 1]])
          }
        }
      }
      if (length(x[[length(x) - 1]]) == 0) { # no Augmenting Path
        return(NA)
      } else { # continue to look for Augmenting Path
        x <- x[-length(x)]
        path <- get_path(x)
        return(path)
      }
    } else { # list ended with Augmenting Path
      for (i in (length(x)/2):1) {
        if (i == length(x)/2) {
          current_edge <- x[[i*2]][c(1,2)]
          path <- c(current_edge, path)
        } else {
          current_edge <- find_by_col(former_edge, x[[i * 2]])
          path <- c(current_edge, path)
        }
        former_edge <- find_by_row(current_edge, x[[i * 2 - 1]])
        path <- c(former_edge, path)
      }
      return(path)
    }
  }
  
  refresh_MM <- function(P_new, MM){ # refresh maximum-matching
    list_P <- list()
    while (length(P_new) > 0) {
      list_P <- c(list_P, list(P_new[c(1, 2)]))
      P_new <- P_new[c(-1, -2)]
    } 
    list_P <- list_P[-1]
    list_MM <- list()
    while (length(MM) > 0) {
      list_MM <- c(list_MM, list(MM[c(1, 2)]))
      MM <- MM[c(-1, -2)]
    }  
    check_P_in_M <- list_P %in% list_MM
    check_M_in_P <- list_MM %in% list_P
    list_MM <- list_MM[!check_M_in_P]
    list_P <- list_P[!check_P_in_M]
    MM <- unlist(c(list_MM, list_P))
    return(MM)
  }
  
  
  ## start calculation ---------
  
  
  MM <- init_match() # maximum-matching, first MM = P
  
  while (TRUE) {
    start_point <- set_entry(MM)
    if (length(start_point) == 0) {
      return(MM)
    } 
    path_list0 <- list(as.vector(start_point))
    path_new <- find_path(start_point, path_list0, MM)
    P_new <- get_path(path_new)
    if (is.na(P_new)[1]) {
      return(MM)
    } else {
      MM <- refresh_MM(P_new, MM)
    }
  }
  
}


## test ---------


# test1 <- matrix(c(1,1,0,0,1,0,1,
#                     1,0,0,0,0,1,0,
#                     0,1,0,0,0,0,0,
#                     0,1,0,1,0,0,0,
#                     0,1,0,1,0,0,0,
#                     0,0,1,1,1,0,0,
#                     0,0,0,0,0,1,0), nrow = 7, ncol = 7, byrow = TRUE)
# test2 <- matrix(c(1,0,1,
#                 0,1,0,
#                 1,0,0),nrow=3, byrow = TRUE)
# > MMHK(test1)
# [1] 3 2 4 4 6 3 7 6 2 1 1 5
# > MMHK(test2)
# [1] 2 2 3 1 1 3

