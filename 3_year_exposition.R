####################################################
# 3 anos de exposição ao redo do censo demográfico #
####################################################

# Estimar a exposição ao longo dos anos civis de 2021-2023, 
# por pessoa que tinha idade inteira A na data do censo de 2022 (1 de agosto de 2022).
# Os resultados são matrizes 100x100 Em e Ef para homens e mulheres, respectivamente.
# Por exemplo:
#   Em %*% [população masculina do censo de 2022 idades 0..99] = 
#   [aprox. exposição masculina de 2021-2022 em 0..99]

rm(list=ls())
library(dplyr)

# Construir funções de sobrevivência l(x) para homens e mulheres, 
# com base no padrão HMD.
# Os multiplicadores calculados não serão muito sensíveis ao padrão escolhido.
# O ponto chave é contabilizar a não sobrevivência em idades avançadas: 
#    uma pessoa de 97,58 anos observada no Censo (t=2022.58) representa mais do 
#    que um ano-pessoa esperado de exposição aos 96 anos durante [2021.0, 2022.0), 
#    por exemplo, porque só vemos os sobreviventes

# log das taxas de mortalidade de um conjunto de países do HMD:
url_data <- "https://raw.githubusercontent.com/mrogon/Topals_mortality_estimates/refs/heads/main/dados/HMDstd.csv"
HMDmort <- read.csv(url_data)

mort <- HMDmort %>%
  mutate( sf = exp(-exp(f)), sm = exp(-exp(m)), # calcula px com base no log(mx)
          lf = c(1, head(cumprod(sf),-1)),# calcula a px acumulada desde o nascimento até cada idade
          lm = c(1, head(cumprod(sm),-1))
  )

lx.male   = approxfun( x=mort$age, y=mort$lm, yleft=1, yright=0)
lx.female = approxfun( x=mort$age, y=mort$lf, yleft=1, yright=0)

# calcula os anos-pessoa de exposição nas idades [A,A+h) x tempos [0,bigT] 
# por indivíduo com exatamente x anos de idade na data do censo bigC
# se 2021.0 é t=0, o Censo de 1º de agosto de 2022 ocorreu em t=1.58 

PY = function(x, A , h=1, bigC=1.58, bigT=3, this.sex='m', dt=.05) {
  tA = bigC - x + A   # data do A-ésimo aniversário
  tgrid = seq(dt/2, bigT-dt/2, dt) # Grade de tempos em que a exposição é possível
  
  this.lx <- lx.male
  if (this.sex == 'f') this.lx <- lx.female
  
  # Para uma pessoa com x anos de idade no tempo bigC, quantos anos-pessoa (p-y) 
  # de exposição nas idades [A, A+n) x tempos [0, T)?
  sum( (tgrid >= tA) * (tgrid < (tA+h)) * this.lx(x-bigC+tgrid) * dt) / this.lx(x)
}

# construir um data frame (no formato long) com todas as combinações de idades
# exatas no censo (em quintos de ano) e faixas etárias de um ano.
# As variáveis expos conterão os anos-pessoa esperados vividos nas idades [A,A+1) 
# ao longo dos anos [0,3], por pessoa de x anos de idade na data do Censo de 1,58 

census.age  = seq(-1.9, 99.9,.20)   # idades possíveis no tempo C
age.group   = 0:99

df = expand.grid( x=census.age, A=age.group)

for (i in 1:nrow(df)) {
  df[i,'m.expos'] = PY(x=df$x[i], A=df$A[i], this.sex='m')
  df[i,'f.expos'] = PY(x=df$x[i], A=df$A[i], this.sex='f')
}  

# Construir uma matriz multiplicadora para converter populações do Censo 
# por idade simples para exposição no período. 
# Agregue idades censitárias exatas em grupos etários inteiros e, então,
# calcule a exposição média. 
# Por exemploEm [100x100] %*% male.census.pop [100x1] = exposição no período masculina [100x1]

tmp = df %>% 
  mutate(intx = floor(x)) %>% 
  group_by(intx,A) %>% 
  summarize(f.expos=mean(f.expos), m.expos=mean(m.expos))

Em  = matrix(round(tmp$m.expos,2), nrow=100, 
             dimnames=list(paste0('A',0:99),paste0('x',-2:99)))


Ef  = matrix(round(tmp$f.expos,2), nrow=100, 
             dimnames=list(paste0('A',0:99),paste0('x',-2:99) ))

# Colapse as primeiras 3 colunas (floor(x) = -2, -1, 0) através da soma. 
# Isso equivale a assumir que as coortes nascidas nos dois anos após o Censo 
# terão o mesmo tamanho que as de 0 anos de idade atuais.


W = rbind( cbind(1, matrix(0, 2, 99)) , 
           diag(100))

Em = Em %*% W
Ef = Ef %*% W

# Verificando os resultados para a população de Sambaiba-MA:

url_data <- "https://raw.githubusercontent.com/mrogon/Topals_mortality_estimates/refs/heads/main/dados/Pop_sambaiba_MA.csv"
sambaiba <- read.csv(url_data, dec=".", sep=';')

N = sambaiba %>% filter(sexo=='Homens') %>% select(pop)
N = N$pop[-101]
N[is.na(N)] <- 0

X = Em %*% N

plot( 0:99, N*3, type='h', lwd=2)
abline(v= seq(0,100,5), lty=2, col='red')
lines(0:99, X, lwd=3, col='seagreen')

save(Em,Ef, file='EmEf2022.RData')

