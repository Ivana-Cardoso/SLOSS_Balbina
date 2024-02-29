################################################################
### Script para selecionar sítios de estudo a partir de um   ###   
###   conjunto pré-definido                                  ###
################################################################
### Autor: Pavel Dodonov - pdodonov@gmail.com                ###
################################################################
### Laboratório de Ecologia Aplicada à Conservação (LEAC)    ###
### PPG em Ecologia e Conservação da Biodiversidade (PPGECB) ###
### Universidade Estadual de Santa Cruz (UESC)               ###
### Ilhéus - BA - Brasil                                     ###
################################################################
### Copyright: CC BY 3.0 US                                  ###
### https://creativecommons.org/licenses/by/3.0/us/          ###
### Cópias, compartilhamentos e modificações são permitidas, ###
###    contanto que o autor original seja mencionado.        ###
################################################################

### Objetivo: a partir de um conjunto pré-definido de sítios 
###    amostrais, definir um subconjunto de sítios que:
### 1) Respeitem uma distância mínima entre eles;
### 2) Maximizem a variação nos valores de uma variável 
###    explanatória.


## Ivana Cardoso:  I modified Pavel's function (with his permission) to randomly select points in each landscape, respecting a minimum distance 

set.seed(13)

select.site <- function(x, var.site, dist.min, coord.X, coord.Y, Nmax, Nsets=10000, Nmin) {
  site.names <- x[,var.site]
  Nsites <- nrow(x)
  sites <- 1:Nsites
  coords <- as.matrix(x[,c(coord.X, coord.Y)])
  dists <- as.matrix(dist(coords))
  dists.allowed <- dists >= dist.min
  sites.rand <- numeric()
  sites.chosen <- list()
  sites.N <- numeric(Nsets)
  
  for(j in 1:Nsets) {
    site.foo <- sample(sites, 1)
    site.rand <- site.foo
    #check distances for this site
    dists.allowed.now <- dists.allowed[site.foo,]
    #remove sites that are not allowed
    sites.new <- sites[dists.allowed.now]
    #start the iterations
    repeat {
      # randomize next site
      if(length(sites.new) > 1) {
        site.foo <- sample (sites.new, 1)
      } else site.foo <- sites.new
      site.rand <- c(site.rand, site.foo)
      # check if enough sites had been sampled
      if(length(site.rand) == Nmax) break
      # check distances for the new iterat
      dists.allowed.now <- dists.allowed[site.rand,]
      dists.allowed.now <- apply(dists.allowed.now,2,all)
      sites.new <- sites[dists.allowed.now]
      if(length(sites.new) == 0) break
    }
    sites.chosen[[j]] <- sort(site.rand)
    sites.N[j] <- length(site.rand)
    if(j %% 500 == 0)  print(c(j,sites.N[j]))
  }
  ### Agora temos conjuntos de sítios com ao menos dist.min entre eles.
  ### Vamos pegar os conjuntos únicos.
  sites.chosen <- unique(sites.chosen)
  ### Vamos remover os que têm poucos sítios
  sites.chosen2 <- sites.chosen[unlist(lapply(sites.chosen, length)) >= Nmin]
  
}

