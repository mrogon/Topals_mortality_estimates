################################################
# Taxas de mortalidade suavizadas com o Topals #
################################################

library(splines)
library(tidyverse)

#---------------------------------
# Carregando as bases de dados
#---------------------------------

url_data <- "https://raw.githubusercontent.com/mrogon/Topals_mortality_estimates/refs/heads/main/dados/Pop_sambaiba_MA.csv"
sambaiba_pop <- read.csv(url_data, dec=".", sep=';')

url_data <- "https://raw.githubusercontent.com/mrogon/Topals_mortality_estimates/refs/heads/main/dados/obitos_sambaiba_2021-2023.csv"
sambaiba_death <- read.csv(url_data, dec=".", sep=',')

url_data <- "https://raw.githubusercontent.com/mrogon/Topals_mortality_estimates/refs/heads/main/dados/Pop_brasil_2022.csv"
brasil_pop <- read.csv(url_data, dec=".", sep=';')

url_data <- "https://raw.githubusercontent.com/mrogon/Topals_mortality_estimates/refs/heads/main/dados/obitos_brasil_2022-2023.csv"
brasil_death <- read.csv(url_data, dec=".", sep=',')

url_raw <- "https://raw.githubusercontent.com/mrogon/Topals_mortality_estimates/refs/heads/main/dados/EmEf2022.RData"
temp <- tempfile(fileext = ".RData")
download.file(url_raw, destfile = temp, mode = "wb")
load(temp); unlink(temp)

#---------------------------
# Parêmetros necessários:
#---------------------------
age = 0:99
B   = bs( 0:99, knots=c(0,1,10,20,40,70), degree=1 )

#------------------------------------------------------
# Calcula o padrão das taxas de mortalidade por idade
#------------------------------------------------------

# Matriz 100x2 com os padrões femininos e masculinos:

br_death <- brasil_death %>%
  filter(idade<100) %>%
  group_by(sexo, idade) %>%
  summarise(obitos=sum(obitos)/3) %>%
  select(sexo, idade, obitos) %>%
  ungroup()

br_pop <- brasil_pop %>%
  filter(idade<100)

br_dados <- left_join(br_pop, br_death, by=c('sexo', 'idade'))

tmp = br_dados %>%
  mutate(lambda = log(obitos/pop))

BR = matrix(tmp$lambda, ncol=2, dimnames=list(0:99, c('Homens','Mulheres')))

## Suavizando o padrão das taxas com uma spline cubica com nós sequênciais:
basis = bs(0:99, knots=seq(0,99,2))
Proj  = basis %*% solve(crossprod(basis)) %*% t(basis)
BRstd = Proj %*% BR

## Verificando os padrões por sexo:

# Homens:
plot(BR[,1])
lines(BRstd[,1])

# Mulheres:
plot(BR[,2])
lines(BRstd[,2])

#-------------------------------------------------
# Função TOPALS para aplicação aos dados:
# soma(óbitos 2021:2023)/exposição (2021:2023)
#-------------------------------------------------

TOPALS_fit = function( N, D, std,
                       max_age        = 99,
                       knot_positions = c(0,1,10,20,40,70), 
                       smoothing_k    = 1,
                       max_iter       = 20,
                       alpha_tol      = .00005,
                       details        = FALSE) {
  
  require(splines)
  
  ## single years of age from 0 to max_age
  age = 0:max_age
  
  ## B is an Ax7 matrix. Each column is a linear B-spline basis function
  B      = splines::bs( age, knots=knot_positions, degree=1 )
  nalpha = ncol(B) 
  
  ## penalized log lik function
  Q = function(alpha) {
    lambda.hat = as.numeric( std + B %*% alpha)
    penalty    = smoothing_k * sum( diff(alpha)^2 )
    return( sum(D * lambda.hat - N * exp(lambda.hat)) - penalty)
  }
  
  ## expected deaths function
  Dhat = function(alpha) {
    lambda.hat = std + B %*% alpha
    return(  as.numeric( N * exp(lambda.hat) ))
  }      
  
  ## S matrix for penalty
  S = matrix(0,nalpha-1,nalpha) 
  diag(S[, 1:(nalpha-1)]) = -1
  diag(S[, 2:(nalpha)  ]) = +1
  SS = crossprod(S)
  
  #------------------------------------------------
  # iteration function: 
  # next alpha vector as a function of current alpha
  #------------------------------------------------
  next_alpha = function(alpha) {
    dhat = Dhat(alpha)
    M = solve ( t(B) %*% diag(dhat) %*% B + 2*smoothing_k *SS)
    v = t(B) %*% (D - dhat) - 2* (smoothing_k * (SS %*% alpha))
    return( alpha + M %*% v)
  }
  
  ## main iteration:     
  a = rep(0, nalpha)
  
  niter = 0
  repeat {
    niter      = niter + 1
    last_param = a
    a          = next_alpha( a )  # update
    change     = a - last_param
    
    converge = all( abs(change) < alpha_tol)
    overrun  = (niter == max_iter)
    
    if (converge | overrun) { break }
    
  } # repeat
  
  if (details | !converge | overrun) {
    if (!converge) print('did not converge')
    if (overrun) print('exceeded maximum number of iterations')
    
    dhat = Dhat(a)
    covar = solve( t(B) %*% diag(dhat) %*% B + 2*smoothing_k *SS)
    
    return( list( alpha    = a, 
                  covar    = covar,
                  Qvalue   = Q(a),
                  converge = converge, 
                  maxiter  = overrun))
  } else return( a) 
  
} # TOPALS_fit

#----------------------------------------------------
# Prepara os dados de Sambaiba para aplicar a função
#----------------------------------------------------

# 3 anos de exposição masculina:
N = sambaiba_pop %>% filter(sexo=='Homens') %>% select(pop)
N = N$pop[-101]
N[is.na(N)] <- 0

Nm = Em %*% N

plot( 0:99, N*3, type='h', lwd=2)
abline(v= seq(0,100,5), lty=2, col='red')
lines(0:99, Nm, lwd=3, col='seagreen')


# 3 anos de exposição feminina:
N = sambaiba_pop %>% filter(sexo=='Mulheres') %>% select(pop)
N = N$pop[-101]
N[is.na(N)] <- 0

Nf = Ef %*% N

plot( 0:99, N*3, type='h', lwd=2)
abline(v= seq(0,100,5), lty=2, col='red')
lines(0:99, Nf, lwd=3, col='seagreen')

#----------------------------------
# Log taxas de mortalidade (mu):
#----------------------------------

#---------
# Homens:
#---------
this.std <- BRstd[,1]

D <- sambaiba_death %>%
  filter(sexo=='Homens') %>%
  group_by(idade) %>%
  summarise(D=sum(obitos)) %>%
  select(D)

D <- D$D
N <- Nm[,1]

mu = D/N

# Roda a função Totals com os argumentos D, N e Std:
fit = TOPALS_fit(N, D, std=this.std, details=TRUE)

# log das taxas suavizadas (logmx):
L = this.std + B %*% fit$a   

# Erro-padrão para o log das taxas suavisadas:
se_logmx = sqrt( diag (B %*% fit$covar %*% t(B)) )

# Intervalo de 95% de confiança (assintótico) para o log das taxas:
Q5.logmx  = L -1.96 * se_logmx
Q95.logmx = L +1.96 * se_logmx

# Vizualizar o resultado no gráfico:
plot(age, log(mu), ylim=c(-10, 1), main="Log(mx) Homens, Sambaiba-MA (2022)")
lines(age, this.std, lwd=3, col='grey')
lines(age, L, lwd=3, col='brown')
lines(age, Q5.logmx, lwd=2, col='brown', lty  = 3)
lines(age, Q95.logmx, lwd=2, col='brown', lty  = 3)


#---------
# Mulheres:
#---------
this.std <- BRstd[,2]

D <- sambaiba_death %>%
  filter(sexo=='Mulheres') %>%
  group_by(idade) %>%
  summarise(D=sum(obitos)) %>%
  select(D)

D <- D$D
N <- Nf[,1]

mu = D/N

# Roda a função Totals com os argumentos D, N e Std:
fit = TOPALS_fit(N, D, std=this.std, details=TRUE)

# log das taxas suavizadas (logmx):
L = this.std + B %*% fit$a   

# Erro-padrão para o log das taxas suavisadas:
se_logmx = sqrt( diag (B %*% fit$covar %*% t(B)) )

# Intervalo de 95% de confiança (assintótico) para o log das taxas:
Q5.logmx  = L -1.96 * se_logmx
Q95.logmx = L +1.96 * se_logmx

# Vizualizar o resultado no gráfico:
plot(age, log(mu), ylim=c(-10, 1), main="Log(mx) Mulheres, Sambaiba-MA (2022)")
lines(age, this.std, lwd=3, col='grey')
lines(age, L, lwd=3, col='brown')
lines(age, Q5.logmx, lwd=2, col='brown', lty  = 3)
lines(age, Q95.logmx, lwd=2, col='brown', lty  = 3)


