library(data.table)
library(ggplot2)


one_sample = function(sample, mu0, alt) {
    mean_x = mean(sample)
    sd_x = sd(sample)
    n_x = length(sample)
    df_x = n_x - 1
    t = sqrt(n_x)*(mean_x-mu0)/sd_x
    if (alt == 0) {
        p_value = pt(t, df_x)
    } else if (alt == 1) {
        p_value = pt(t, df_x, lower.tail = FALSE)        
    } else {
        p_value = 2*pt(abs(t), df_x, lower.tail = FALSE) 
    }
    list(stat = t, p = p_value)
}

independent = function(sample, sample2, alt) {    
    mean_x = mean(sample)
    var_x = var(sample)
    n_x = length(sample)
    mean_y = mean(sample2)
    var_y = var(sample2)
    n_y = length(sample2)
    df = n_x + n_y - 2
    t = (mean_x-mean_y)/sqrt((n_x-1)*var_x+(n_y-1)*var_y)*sqrt((n_x*n_y*df) / (n_x+n_y))
    if (alt == 0) {
        p_value = pt(t, df)
    } else if (alt == 1) {
        p_value = pt(t, df, lower.tail = FALSE)        
    } else {
        p_value = 2*pt(abs(t), df, lower.tail = FALSE) 
    }
    list(stat = t, p = p_value)
}

dependent = function(sample, sample2, alt) {
    d = sample - sample2
    mean_d = mean(d)
    sd_d = sd(d)
    n_d = length(d)
    df_d = n_d - 1
    t = sqrt(n_d)*mean_d/sd_d
    if (alt == 0) {
        p_value = pt(t, df_d)
    } else if (alt == 1) {
        p_value = pt(t, df_d, lower.tail = FALSE)        
    } else {
        p_value = 2*pt(abs(t), df_d, lower.tail = FALSE) 
    }
    list(stat = t, p = p_value)
}

test_t_student <- function(test,data) {
    if (test == 1) {
    column = readline(prompt = "Give the column name:\n ")
    mean = readline(prompt = "Give the theoretical mean:\n ")
    alt <- as.numeric(readline(prompt = "Choose the alternative hypothesis (0 - less, 1 - greater, 2 - two sided):\n "))     
    wilk = shapiro.test(as.numeric(data[[column]]))

    print(wilk)    
    one_sample(data[[column]],as.numeric(mean),alt)    

  } else if (test == 2){
    column_1 = readline(prompt = "Give the first column name:\n ")
    column_2 = readline(prompt = "Give the second column name:\n ") 
    alt <- as.numeric(readline(prompt = "Choose the alternative hypothesis (0 - less, 1 - greater, 2 - two sided):\n "))

    wilk1 = shapiro.test(as.numeric(data[[column_1]]))
    wilk2 = shapiro.test(as.numeric(data[[column_2]])) 
    variance = var.test(data[[column_1]], data[[column_2]])

    print(wilk1)
    print(wilk2)
    print(variance)
    independent(data[[column_1]],dane[[column_2]],alt)

  } else if (test == 3){
    column_1 = readline(prompt = "Give the first column name:\n ")
    column_2 = readline(prompt = "Give the second column name:\n ") 
    alt <- as.numeric(readline(prompt = "Choose the alternative hypothesis (0 - less, 1 - greater, 2 - two sided):\n "))


    wilk1 = shapiro.test(as.numeric(data[[column_1]]))
    wilk2 = shapiro.test(as.numeric(data[[column_2]])) 
    variance = var.test(data[[column_1]], data[[column_2]])

    print(wilk1)
    print(wilk2)
    print(variance)
    independent(data[[column_1]],dane[[column_2]],alt)  

  } else {
   print("Incorrect value. Choose 1, 2 or 3.\n")
  }
}


lm2 = function(y, int, variable_names, data1) {
      if (int == 0) {
        X <- as.matrix(do.call(cbind, lapply(variable_names, function(col) data1[[col]])))
      } else {
        X <- as.matrix(cbind(rep(1, length(y)), do.call(cbind, lapply(variable_names, function(col) data1[[col]]))))
      }
      y <- as.matrix(y)

      beta <- solve(t(X) %*% X) %*% t(X) %*% y

      if (int == 0) {
        res <- y - X %*% beta
      } else {
        res <- y - X[, 1] * beta[1] - X[, -1, drop = FALSE] %*% beta[-1]
      }

      S2e <- t(res) %*% res / (nrow(y) - ncol(X))
      Se <- sqrt(S2e)

      seBeta <- sqrt(diag(c(S2e) * solve(t(X) %*% X)))
      t <- beta / seBeta

      p.value.t <- 2 * pt(abs(t), nrow(y) - ncol(X), lower.tail = FALSE)

      if (int == 0) {
        R2 <- 1 - t(res) %*% res / sum((y^2))
      } else {
        R2 <- 1 - t(res) %*% res / sum((y - mean(y))^2)
      }

      if (int == 0) {
        F <- R2 / (1 - R2) * (nrow(y) - ncol(X)) / ncol(X)
        p.value.f <- pf(F, ncol(X), nrow(y) - ncol(X), lower.tail = FALSE)
      } else {
        F <- R2 / (1 - R2) * (nrow(y) - ncol(X)) / (ncol(X) - 1)
        p.value.f <- pf(F, ncol(X) - 1, nrow(y) - ncol(X), lower.tail = FALSE)
      }

      list(est = beta, p1 = p.value.t, r2 = R2, f = F, p2 = p.value.f)
    }

anova = function(y, alpha) {
  k = sort(unique(alpha))

  n.i = numeric()
  y.i = numeric()
  for (i in k) {
    n.i[i] = length(alpha[alpha == i])
    y.i[i] = mean(y[alpha == i])
  }

  m.y = mean(y)

  SSb = sum(n.i * (y.i - m.y)^2)
  S2b = SSb / (length(k) - 1)

  SSt = sum((y - m.y)^2)
  S2t = SSt / (length(y) - 1)

  SSw = SSt - SSb
  S2w = SSw / (length(y) - length(k))

  F = S2b / S2w
  p = pf(F, length(k) - 1, length(y) - length(k), lower.tail = FALSE)



  differences = NULL
  for (i in 1:(length(k)-1)) {
    for (j in (i+1):length(k)) {
      diff = abs(y.i[i] - y.i[j])
      if (diff > 1.96 * sqrt(S2w / n.i[i] + S2w / n.i[j])) {
        differences = c(differences, paste("Groups", k[i], "i", k[j]))
      }
    }
  }

  list(statistics = F, p.value = p, differences = differences)

}



select_test <- readline("What test do you need to perform (t-test, regression, anova)?:\n ")


if (tolower(select_test) == "t-test") {
    path = readline(prompt = "Give the path to file:\n ")
    data = fread(path)
    test = readline(prompt = "What kind of t-test you need to perform:
    1 - for one sample
    2 - for two independt samples
    3 - for two dependt samples\n")
    test_t_student(test,data)

} else if (tolower(select_test) == "regression") {
    path1 = readline(prompt = "Give the path to file:\n ")
    data1 = fread(path1)
    head(data1)

    int = readline(prompt = "Do you need an intercept?
    1 - yes
    0 - no\n")

    interpreted = readline(prompt = "Give the name of the interpreted column:\n")
    variable_names <- list()
    variable_number <- as.integer(readline("Give the number of interpreted variables:\n "))
    for (i in 1:variable_number) {
      variable_name <- readline(paste("Give the name of the explanatory column", i, ":\n "))
      variable_names[[i]] <- variable_name
    }
    y=dane1[[interpreted]] 


    plots <- list()
    for (i in 1:length(variable_names)) {
      variable_name <- variable_names[[i]]
      plot <- ggplot(data1, aes_string(x = variable_name, y = interpreted)) +
        geom_point() +
        labs(x = variable_name, y = "y") +
        ggtitle(paste("Relationship between", variable_name, "and interpreted variable"))
      plots[[i]] <- plot
    }
    for (i in 1:length(plots)) {
      print(plots[[i]])
    }

    
    if (is.numeric(y)) {
      print(" y are numeric variables")
    } else {
      print("y doesn't have numeric variables")
    }
    for (i in 1:length(variable_names)) {
      variable_name <- variable_names[i]]
      if (is.numeric(data1[[variable_name]])) {
        print(paste("Variables in column", variable_name, "are numeric."))
      } else {
        print(paste("Variables in column", variable_name, "are not numeric"))
      }
    }
    
    
    column_length <- sapply(variable_names, function(col) length(data1[[col]]))
    equal  <- all(column_length == length(y))
    if (equal) {
      print("The number of elements in the columns is equal to the number of elements in the variable y.")
    } else {
      print("The number of elements in the columns is differ from the number of elements in the variable y.")
    }


    lm2(y,int,variable_names,data1)
    
    
} else if (tolower(select_test) == "anova") {
    path2 = readline(prompt = "Give the path to file:\n ")
    data2 = fread(path2)
    head(data2)
    y = readline(prompt = "Give the name of the interpreted column\n")
    alpha = readline(prompt = "Give the name of the column with the factor\n")
    y = data2[[y]]
    alpha = data2[[alpha]]

    anova(y, alpha)
    
} else {
  cat("Invalid input.\n")
}