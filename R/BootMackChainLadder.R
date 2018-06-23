#' @import ChainLadder
#' @import tidyverse
#' @import magrittr
#' @import dplyr
#' @import purrr
#' @import mvtnorm
#' @import ggplot2
#' @import gridExtra
#' @import dplyr
NULL


##### Basic functions                                                                #####
.diag        <- function(Triangle){
# sympli returns the diagonal of the Triangle
  unname(rev(diag(Triangle[nrow(Triangle):1,])))
}
.DF          <- function(Triangle){
# Simply returns chain-ladder simples developpement factors.
  n = dim(Triangle)[1]
  Triangle2 <- Triangle
  Triangle2[ col(Triangle2) == n - row(Triangle2)  + sum(!is.na(Triangle2[,n]))  ] <- NA
  return(unname(c(colSums(Triangle[,2:n],na.rm=TRUE)/colSums(Triangle2[,1:(n-1)],na.rm=TRUE),1)))
}
.ultimates   <- function(Triangle,df = .DF(Triangle)){
  # simply returns chain-ladder ultimates.
  .diag(Triangle) * cumprod(rev(df))
}
.ibnr        <- function(Triangle,df = .DF(Triangle)){
  # simply returns the IBNR from classical chain-ladder model
  .diag(Triangle) * (cumprod(rev(df))-1)
}
.DFIndiv     <- function(Triangle){
  # Simply returns the Triangle of individual developpements factors.
  unname((cbind(Triangle,NA)/cbind(NA,Triangle))[,2:(dim(Triangle)[[1]]+1)])
}
.sigma       <- function(Triangle,df=.DF(Triangle),DFIndiv = .DFIndiv(Triangle)){
  # Returns Mack's estimate for the variances parameters sigma_j
  # Returns a line vector, indexed by Triangle columns.
  n = dim(Triangle)[[1]]
  ecart = DFIndiv - matrix(df,nrow=n,ncol=n,byrow=TRUE)
  rez = Triangle*ecart^2
  sigma = colSums(rez,na.rm=TRUE) / ((n-2):(-1))

  # the last value is approximated by mack's sugestion :
  sigma = c(sigma[1:(n-2)],
            min(sigma[n-2]^2/sigma[n-3],
                sigma[n-2],
                sigma[n-3])
  )
  if(is.nan(sigma[n-1])){
    sigma[n-1] <- 0
  }
  # we return not sigma^2 but just sigma
  return(unname(sqrt(sigma)))
}
.rho         <- function(triangles,sigma = lapply(triangles,.sigma),indivF = lapply(triangles,.DFIndiv),f = lapply(triangles,.DF)){
  # estimateur de la matrice de correlation par colonne entre les differents triangles :
  # retourne un array [N - 1 , NbTriangles,Nbtriangles] contenant les matrices de correlations entre les colonnes des triangles d'input.
  # on attend dans triangles une liste de triangles

  # construction de l'objet array que l'on va retourner :
  I = dim(triangles[[1]])[1]
  rho <- array(NA,dim = c(I-1,length(triangles),length(triangles)))
  omega <- rho

  # # Precomputage des differentes quantitee d'interet pour nos triangles :
  # sigma <- lapply(triangles,.sigma)
  # indivF <- lapply(triangles,.DFPassageIndiv)
  # f <- lapply(triangles,.DFPassage)

  # Construison rho :
  for(n in seq_along(triangles)){
    for(m in seq_along(triangles)){
      if(n==m) {
        # si n = m, on est sur le meme triangle
        rho[,n,m] = 1

      } else {
        # si n != m, alors on va construire une matrice de correlation :

        for(j in 1 : (I - 2)){
          # On commencer par contruire w_j^2
          omegaup = 0
          omegadown = 0
          for(i in 1:(I-(j-1)-1)){
            # Partie hautte e la fraction
            omegaup <- omegaup + sqrt(triangles[[n]][i,j] * triangles[[m]][i,j])

            # partie basse de la fraction :
            somme = 0
            for(k in 1:(I-(j-1)-1)){
              somme = somme + triangles[[m]][k,j]
            }
            omegadown <- omegadown + triangles[[n]][i,j] * somme
          }
          omegaj <- omegaup^2 / omegadown

          # Puis on construit c_j
          cj <- 1/(I -(j-1)-2+omegaj)

          # Enfin on construit la partie somme de rho_j_n_m :
          somme = 0
          for(l in 1:(I-(j-1)-1)){
            somme <- somme + sqrt(triangles[[n]][l,j] * triangles[[m]][l,j])*(indivF[[n]][l,j] - f[[n]][j])*(indivF[[m]][l,j] - f[[m]][j])/(sigma[[n]][j]*sigma[[m]][j])
          }
          # puis on agrege ces differents resultats :
          rho[j,n,m] <- cj * somme

        }

        # Il nous reste a estimer le dernier facteur par une prolongation polynomiale (proposee par MAK 1997, comme poru la vairance : ):
        rho[I-1,n,m] <- min( (rho[I-2,n,m]*sigma[[n]][I-2]*sigma[[m]][I-2])^2/(rho[I-3,n,m]*sigma[[n]][I-3]*sigma[[m]][I-3]),
                             abs(rho[I-3,n,m]*sigma[[n]][I-3]*sigma[[m]][I-3]),
                             abs(rho[I-2,n,m]*sigma[[n]][I-2]*sigma[[m]][I-2]))

      }}}
  # On retourne un array (j,n,m) avec j l'indice de la colonne et n et m les indices des triangles.
  # On en profite juste pour metrte des zero la ou il a pas pu les calculer :
  rho[is.nan(rho)] <- 0
  return(rho)
}
.residuals   <- function(Triangle,centered = FALSE,DF=.DF(Triangle),DFIndiv=.DFIndiv(Triangle),sigma=.sigma(Triangle)){
  # estimation of mack's residuals, based on England & Verral 2007
  # returns a Triangle of same size as the input, but with NA on the diagonal (no residuals for the diagonal)
  residus <- Triangle
  n <- dim(Triangle)[1]

  for(i in 1:n){ for(j in 1:n){
    residus[i,j] <- sqrt((n-j)/(n-j-1))*sqrt(Triangle[i,j]) * (DFIndiv[i,j] - DF[j]) / sigma[j]
  }}
  residus[is.nan(residus)] <- 0
  if(centered){
    residus <- apply(residus,2,function(x){return(x -mean(x,na.rm=TRUE))})
  }
  return(residus)
}
.sampling    <- function(Triangle,N=100,seuil = NA,zonnage=FALSE){
  # Fonction qui bootstrap un Triangle de positions
  # On fabrique une liste de Triangle de r?sidus bootstrapp?s :

  n = length(Triangle)
  I = dim(Triangle)[1]

  # Gestion du seuil d'exclusion des residus :
  Triangle2 <- Triangle
  if(!is.na(seuil)){
    Triangle2[Triangle[(!is.na(Triangle))]>seuil] <- NA
  }
  if(zonnage){
    # Il faut alors zonner les residus. les morceaux de zonnage sont fixes et ne sont pas modifiables
    morceaux = list(c(1,2),seq(3,9),seq(10,12),seq(13,I))

    samples <- lapply(1:length(morceaux),function(i){
      prob = !is.na(Triangle2[,morceaux[[i]]])
      positions <- matrix(1:n,nrow=I)
      return(lapply(1:N,function(x){
        matrix(sample(positions[,morceaux[[i]]],
                      replace = TRUE,
                      prob = prob,
                      size = length(prob))
               ,nrow=I)
        }))
    })

    rez <- lapply(1:N,function(i){
      return(do.call(cbind,lapply(samples,function(x){x[[i]]})))
    })

  } else {
    prob = !is.na(Triangle2)
    positions <- matrix(1:n,nrow=I)
    rez <- sample(positions,replace = TRUE,prob = prob,size = N * n)
    rez <- lapply(1:N,function(.x) {matrix(rez[(n*(.x-1)+1):(n*.x)],nrow=I)})
  }

  rez <- lapply(rez,function(x) {
    x[is.na(Triangle)] <- NA
    return(x)
  })
  return(rez)
}
.applysample <- function(sample,Triangle){
  # Fonction d'application de cet sampling :
  return(ChainLadder::as.triangle(
    matrix(as.vector(Triangle)[as.vector(sample)],
           nrow=dim(Triangle))))
}
.DFPond      <- function(DFIndiv,ponderation){
  n <- dim(ponderation)[1]
  ponderation[col(ponderation)==n-row(ponderation)+1] <- NA
  rez = colSums(ponderation*DFIndiv,na.rm=TRUE)/colSums(ponderation,na.rm=TRUE)
  rez[n] <- 1
  return(rez)

}
.newDFPond   <- function(NyDFIndiv,DFIndiv,Triangle){
  n <- dim(Triangle)[1]
  DFIndiv[col(DFIndiv) == n - row(DFIndiv) + 1] <- c(NyDFIndiv,1)
  rez = colSums(Triangle*DFIndiv,na.rm=TRUE)/colSums(Triangle,na.rm=TRUE)
  rez[n] <- 1
  return(rez)
}
.formatOutput <- function(n,Triangle,diag,DF,sigma,residuals,DFBoot,Ultimates,IBNR,NyCum,NyInc,NyDF,NyUltimates,NyIBNR){
  # output formatting :
  NyCum       <- do.call(rbind,lapply(NyCum,rev))
  NyInc       <- do.call(rbind,lapply(NyInc,rev))
  NyUltimates <- do.call(rbind,lapply(NyUltimates,rev))
  NyIBNR      <- do.call(rbind,lapply(NyIBNR,rev))
  NyDF        <- do.call(rbind,NyDF)
  DFBoot      <- do.call(rbind,DFBoot)

  names(Ultimates)      <- rownames(Triangle)
  names(IBNR)           <- rownames(Triangle)
  names(diag)           <- rownames(Triangle)
  colnames(NyCum)       <- rownames(Triangle)[2:(n)]
  colnames(NyInc)       <- rownames(Triangle)[2:(n)]
  colnames(NyUltimates) <- rownames(Triangle)
  colnames(NyIBNR)      <- rownames(Triangle)
  colnames(NyDF)        <- colnames(Triangle)
  colnames(DFBoot)      <- colnames(Triangle)

  result <- list(
    Latest = rev(diag),
    DF = DF,
    sigma = sigma,
    residuals = residuals,
    DFBoot = DFBoot,
    Ultimates = Ultimates,
    IBNR = IBNR,
    NyCum  = NyCum,
    NyInc =NyInc,
    NyDF = NyDF,
    NyUltimates = NyUltimates,
    NyIBNR = NyIBNR
  )
  class(result) <- c("BootMackChainLadder",class(result))
  return(result)
}
.passage.a.lultime      <- function(tri.ds,cad){
  # passage a l'ultime :
  I <- dim(tri.ds)[[1]]
  tri.ds.ult <- tri.ds
  for(i in 1:I){
    for (j in 1:(I-i+1)){
      tri.ds.ult[i,j] <- tri.ds[i,j]/cad[i+j-1]
    }
  }

  # On retourne le triangle :
  return(tri.ds.ult)
}
.erreur.passage.ultime  <- function(tri.ds.ult,tri.si,tri.ds) {
  # Clacul de l'erreur du au passage a l'ultime :
  var1 <- tri.ds.ult %>% .diag %>% sum
  var2 <- tri.ds %>% sum(na.rm=TRUE)
  var3 <- tri.si %>% .ibnr %>% sum

  return(var1 - var2 - var3)
}
.check.triangle         <- function(triangle,dim=dim(triangle)){
  # FOnction de verification des dimentions d'un triangle :
  lignes <- row.names(triangle) %>%
    as.numeric %>%
    {add(.,-.[1])} %>%
    {prod(. == (0:(length(.)-1)))} # Retourne TRUE si les lignes sont bonnes, faux sinon.

  collonnes <- colnames(triangle)  %>%
    as.numeric %>%
    {add(.,-.[1])} %>%
    {prod(. == (0:(length(.)-1)))}  # Retourne TRUE si les collonnes sont bonnes, faux sinon.

  return(lignes & collonnes) # si le retour est faux, il faut fixer le triangle.
}
.fix.triangle.dimention <- function(triangle,taille=27,replaceNAbyZero = FALSE,deja_cumule=TRUE){



  # des zero a la place des NA :
  triangle[is.na(triangle)] <- 0

  if(deja_cumule){
    # Passage en incremental si et seulement si ca renvois pas une erreur :
    lastfunc <- function(x){
      incr2cum(x,na.rm=TRUE)
    }
    triangle <- tryCatch(expr={cum2incr(triangle)},error=function(e){
      lastfunc <<- function(x){x}
      triangle
    })
  }




  tri <- matrix(NA,nrow(triangle),taille)
  tri2 <- matrix(NA,taille,taille)

  # Occupons nous d'abords de completer les colonnes :
  for(i in 1:taille){
    if(as.character(i-1) %in% colnames(triangle)){
      tri[,i] <- triangle[,colnames(triangle) == as.character(i-1)]
    } else{
      tri[,i] <- rep(NA,nrow(triangle))
    }
  }

  # puis completons les lignes :
  zero = 1991
  for (i in 1:taille ){
    if((i-1) %in% (as.numeric(row.names(triangle))-zero)){
      tri2[i,] <- tri[which((as.numeric(row.names(triangle))-zero) %in% (i-1)),]
    } else {
      tri2[i,] <- rep(NA,taille)
    }
  }

  # on reconstruit le triangle, on remet les noms en place et on retourne :
  tri <- as.triangle(tri2)

  # On remplace tout les NA sous la diagonale par des 0 :
  for(i in 1:taille){
    for(j in 1:(taille-i+1)){
      if(is.na(tri[i,j])){
        tri[i,j] = 0
      }
    }
  }
  for(i in 1:taille){
    for(j in 1:taille){
      if(i+j>taille+1){
        tri[i,j] <- NA
      }
    }
  }

  # on remet en place les noms de lignes et de colonnes :
  row.names(tri) <- as.character(1991:(1991+taille-1))
  colnames(tri) <- as.character(0:(taille-1))


  if(deja_cumule){
    tri <- lastfunc(tri)
  }
  # on cumule et on renvois :
  return(tri)
}
.cumsd <- function(x){
  return(sapply(seq_along(x), function(k, z) sd(z[1:k]), z = x))
}
.cumcorel <- function(x){
  noms = names(x)
  nomMat = map(1:ncol(x),~map(1:ncol(x),function(y){paste0(noms[.x],"-",noms[y])}))
  cor = list()
  for(i in 1:nrow(x)){
    cor[[i]] = cor(x[1:i,]) %>% as.vector
  }
  cor %>%
  {do.call(rbind.data.frame,.)} %>%
    set_names(nomMat %>% unlist) %>%
    return
}
##### Boot Mack Chain Ladder                                                         #####
#' Boot Mack Chain Ladder model
#'
#' This function implement a simple bootstrap of the residuals from the mack model with a one-year reserving risk point of view.
#'
#' @param Triangle A simple triangle from che Chain-ladder package.
#' @param B numeric. The number of bootstrap samples you want
#' @param distNy character. Distribution of next-year incremental payments. Either "normal" (default) or "residuals"
#' @param seuil numeric. Seuil for residuals exclusion. A value of NA (default) will prevent
#' @param zonnage logical. Do you want to zonne the residuals ?
#'
#' @return A BootMackChainLadder object with a lot of information about the bootstrapping. You can plot it, print it, str it to extract information.
#' @export
#'
#' @import magrittr
#'
#' @examples
#' data(ABC)
#' BootMackChainLader(Triangle = ABC, B = 100, distNy = "residuals", seuil = 2)
BootMackChainLadder <- function(Triangle,B=100,distNy="normal",seuil=NA,zonnage=FALSE){

  if(!(distNy %in% c("normal","residuals"))){stop("DistNy Parameter must be 'normal' (classical MW) or 'residuals'")}
  # Cf "One-year reserve risk including a tail factor : closed formula and bootstrapp aproach, 2011"

  # First step : Mack model.
  n         <- dim(Triangle)[1]
  diag      <- rev(.diag(Triangle))
  DF        <- .DF(Triangle)
  DFIndiv   <- .DFIndiv(Triangle)
  sigma     <- .sigma(Triangle,DF,DFIndiv)
  residuals <- .residuals(Triangle,centered = TRUE,DF,DFIndiv,sigma)
  Ultimates <- .ultimates(Triangle,DF)
  IBNR      <- .ibnr(Triangle,DF)

  # Step 2 : Resampling residuals.
  samples          <- .sampling(residuals,B,seuil = seuil,zonnage = zonnage)
  sampledResiduals <-lapply(samples,.applysample,residuals)
  # Step3 : Calculation of booted estimators
  DFIndivBoot      <- lapply(sampledResiduals,function(.x){t(t(.x*sqrt(t(c(sigma^2,0)/t(Triangle))))+DF)})
  DFBoot           <- lapply(DFIndivBoot,.DFPond,Triangle) # ponderated by the original triangle !

  # Step 5 : Simulation of NY payments by a standard normal OR by the residuals... taking into acount process error
  if(distNy=="normal"){
    NyCum <- lapply(DFBoot,function(.x){rnorm(n = n-1,mean = (diag * .x)[1:(n-1)],sd =   sqrt((diag * c(sigma^2,0))[1:(n-1)]))})
  }
  if(distNy=="residuals"){
    # To simulate the same thing with residuals assuming they are all i.i.d between developpements years, we could do :
    if(!is.na(seuil)){ # FIrst, let's discard non-selected residuals.
      residuals[residuals[(!is.na(residuals))]>seuil] <- NA
    }
    if(!zonnage){
      samples <- sample(residuals[!is.na(residuals)],size=(n-1)*B,replace = TRUE)
      samples <- lapply(1:B,function(.x) {samples[((n-1)*(.x-1)+1):((n-1)*.x)]})
      NyCum <- mapply(function(.y,.x){(.y *sqrt((diag * c(sigma^2,0))[1:(n-1)]))+(diag * .x)[1:(n-1)]},samples,DFBoot,SIMPLIFY = FALSE)
    } else {
      # If zonnig of residuals is needed, the following zonning will be implemented :
      # first from collumn 1 to column 2
      # then from column 3 to column column 9
      # then from 10 to 12
      # then 13+
      # not implemented yet !
      morceaux = list(c(1,2),seq(3,9),seq(10,12),seq(13,n-1))
      samples = list()

      samples <- lapply(1:length(morceaux),function(i){
        x <- sample(residuals[,morceaux[[i]]][!is.na(residuals[,morceaux[[i]]])],size=length(morceaux[[i]])*B,replace=TRUE)
        x <- lapply(1:B,function(.x) {x[(length(morceaux[[i]])*(.x-1)+1):(length(morceaux[[i]])*.x)]})
        return(x)
      })


      # for(i in 1:length(morceaux)){
      #   samples[[i]] <- sample(residuals[,morceaux[[i]]][!is.na(residuals[,morceaux[[i]]])],size=length(morceaux[[i]])*B,replace=TRUE)
      #   samples[[i]] <- lapply(1:B,function(.x) {samples[[i]][(length(morceaux[[i]])*(.x-1)+1):(length(morceaux[[i]])*.x)]})
      # }


      samples2 <- lapply(1:B,function(i){
        plop <- samples[[1]][[i]]
        for(j in 2:length(morceaux)){
          plop <- c(plop,samples[[j]][[i]])
        }
        return(plop)
      })


      # samples2 <- list()
      # for(i in 1:B){
      #   plop <- samples[[1]][[i]]
      #   for(j in 2:length(morceaux)){
      #     plop <- c(plop,samples[[j]][[i]])
      #   }
      #   samples2[[i]] <- plop
      # }

      NyCum <- mapply(function(.y,.x){(.y *sqrt((diag * c(sigma^2,0))[1:(n-1)]))+(diag * .x)[1:(n-1)]},samples2,DFBoot,SIMPLIFY = FALSE)


    }

  }

  # NyCum is list of B vectors of size n-1, containing the next year cumulative payments, from origin n to 1 ( decreasing order)

  #Step 6 : Calculation of Ny estimators
  NyInc            <- lapply(NyCum,function(.x){.x - diag[1:(n-1)]}) # Coresponding incremnts
  NyDFIndiv        <- lapply(NyCum,function(.x){.x/diag[1:(n-1)]}) # Coresponding individual DF's
  NyDF             <- lapply(NyDFIndiv,.newDFPond,DFIndiv,Triangle) # coredponing DF
  NyUltimates      <- mapply(function(.x,.y){c(.x*rev(cumprod(rev(.y[2:n]))),Triangle[1,n])},NyCum,NyDF,SIMPLIFY = FALSE)
  NyIBNR <- lapply(NyUltimates,function(.x){.x-diag})

  rez <- .formatOutput(n,Triangle,diag,DF,sigma,residuals,DFBoot,Ultimates,IBNR,NyCum,NyInc,NyDF,NyUltimates,NyIBNR)
  return(rez)
}

#' mean.BootMackChainLadder
#'
#' Calculate mean statiscics from the bootstraped mack model
#'
#' @param x A BootMackChainLadder Object
#'
#' @return Three data.frames containing data per origin year ("ByOrigin"), by developpement year ("ByDev) and globaly ("Totals)
#' @export
#'
#' @examples
#' data(ABC)
#' BMCL <- BootMackChainLader(Triangle = ABC, B = 100, distNy = "residuals", seuil = 2)
#' mean(BMCL)
mean.BootMackChainLadder    <- function(x,...){

  # the purpose of this function is to return means of everything BMCL returned.
  ByOrigin = data.frame(
    Latest = x$Latest,
    Ultimates = x$Ultimates,
    IBNR = x$IBNR,
    NyCum  = c(x$Ultimates[1],colMeans(x$NyCum)),
    NyInc =c(0,colMeans(x$NyInc)),
    NyUltimates = colMeans(x$NyUltimates),
    NyIBNR = colMeans(x$NyIBNR)
  )
  row.names(ByOrigin) = rev(row.names(ByOrigin))

  ByDev = data.frame(
    DFBoot = colMeans(x$DFBoot),
    NyDF = colMeans(x$NyDF)
  )

  Totals = colSums(ByOrigin)


  return(list(ByOrigin = ByOrigin, ByDev = ByDev,Totals = Totals))
}
#' CDR.BootMackChainLadder
#'
#' Calculate the one-year Claim developpement results from the bootstrapped mack model.
#'
#' @param BMCL A BootMackChainLadder Object
#'
#' @return A data.Frame with one-year results from the bootstrap.
#' @export
#'
#' @examples
#' data(ABC)
#' BMCL <- BootMackChainLader(Triangle = ABC, B = 100, distNy = "residuals", seuil = 2)
#' CDR(BMCL)
CDR.BootMackChainLadder     <- function(BMCL){
  Totals = data.frame(IBNR = sum(BMCL$IBNR),`CDR(1)S.E.` = sd(rowSums(BMCL$NyIBNR)))
  row.names(Totals) <- "Totals"
  names(Totals) <- c("IBNR","CDR(1)S.E.")
  rbind(setNames(data.frame(IBNR = BMCL$IBNR,`CDR(1)S.E.` = apply(BMCL$NyIBNR,2,sd)),names(Totals)), Totals)

}
#' summary.BootMackChainLadder
#'
#'  Give summary statistics about the bootstrap.
#'
#' @param object A BootMackChainLadder Object
#'
#' @return NULL
#' @export
#'
#' @details
#' The function only print information about the model.
#'
#' @examples
#' data(ABC)
#' BMCL <- BootMackChainLader(Triangle = ABC, B = 100, distNy = "residuals", seuil = 2)
#' summary(BMCL)
summary.BootMackChainLadder <- function(object,...){
  mean <- mean(object)
  cat("This is a BootMackChainLadder model \n\n")
  cat("Detailed results by Origin years : \n")
  print(format(mean$ByOrigin,format="i", nsmall=0, big.mark=","))
  #print(mean$ByDev)
  cat("\n Totals across origin years : \n")
  print(format(t(mean$Totals),format="i", nsmall=0, big.mark=","))
  return(NULL)
}
#' print.BootMackChainLadder
#'
#' @param BMCL
#'
#' @return the BMCL object
#' @export
#'
#' @examples
#' data(ABC)
#' BMCL <- BootMackChainLader(Triangle = ABC, B = 100, distNy = "residuals", seuil = 2)
#' print(BMCL)
print.BootMackChainLadder   <- function(x,...){
  print(summary(x))
  return(NULL)
}


##### Multi Boot Back Chain Ladder                                                   #####
#' MultiBootMackChainLadder
#'
#' The multi boot mack chain ladder algorythme computes todays and next-year common estimates on a portefolio of several triangles, following closely a synchronised version of BootMackChainLadder.
#'
#' @param triangles A List of Triangles objects of the same size.
#' @param B The numebr of boostrap replicates
#' @param distNy The process distribution for next year increments. Either "normal" (default), "residuals.bycolumn","residuals.global" or "residuals". See details.
#' @param names enventual names of the different triangles. IF set to NULL, the procedure will try to get names from the triangles list.
#' @param seuil Eventual exclusions limit for residuals. Set to NA (default) to avoid excluding anything.
#'
#' @details
#'
#' This model uses the fact that the Mack model can be seen as a quasi-glm to found nice residuals. Bootstrap thoose residuals on the upper-left triangle allows to get bootstrap distribution of today's estimatins ( reserves, ultimates, ...).
#'
#' In each bootstrap resample, the function use the specified process distrubiton to simulate next-year payments. If a normal law is used, it follows boumezoued et all and converges to the Mer-wuthrich fromula in the Braun model. If set to residuals, the convergence is there but not to the same result since no specific law is supposed for residuals.
#'
#' @return a MBMCL object containing a list of BMCL objects and a little more.
#' @export
#'
#' @seealso BootMackChainLadder
#'
#' @examples
#' data(ABC)
#' triangles <- list(tri1 = ABC, tri2 = ABC, tri3 = ABC)
#' MultiBootmackChainLadder(triangles,100)
MultiBootMackChainLadder    <- function(triangles,B=100,distNy = "normal", names=NULL,seuil=NA){
  if(!(distNy %in% c("normal","residuals.bycolumn","residuals.global","residuals"))){stop("DistNy Parameter must be 'normal' (classical MW) or 'residuals'")}
  # Cf "One-year reserve risk including a tail factor : closed formula and bootstrapp aproach, 2011"

  # First step : deterministic statistics about the inputed triangles
  n         <- dim(triangles[[1]])[1] # dimention des triangles
  N         <- length(triangles) # nombre de triangles
  diag      <- lapply(triangles,function(.x){rev(.diag(.x))})
  DF      <- lapply(triangles,.DF)
  DFIndiv <- lapply(triangles,.DFIndiv)
  sigma     <- mapply(.sigma,triangles,DF,DFIndiv,SIMPLIFY = FALSE)
  residuals <- mapply(.residuals,triangles,center=TRUE,DF,DFIndiv,sigma,SIMPLIFY = FALSE)
  ibnr      <- mapply(.ibnr,triangles,DF,SIMPLIFY = FALSE)
  rho       <- .rho(triangles,sigma,DFIndiv,DF)
  Ultimates <- mapply(.ultimates,triangles,DF,SIMPLIFY=FALSE)

  # Second step : resampling and calculation of basic statistics on resampled triangles.


  # first step : tweaking residuals to exclude somes.
  # the exclusions needs to be done sychronised. So our proposal is to replace positions that are excluded by positions from the same original column.

    # Let's calculate excluded positions :

  if(!is.na(seuil)){
    # positions that pose problem :
    excluded_positions <- lapply(residuals,function(x){which(abs(x)>seuil)})

    # on thoose positions, we have to replace the residuals that pose problem by one another from the same collumn.
    residuals <- mapply(function(res,pos){
      for(i in 1:length(pos)){
        col = unique(floor(pos[i]/n)+1)
        line <- which(abs(res[((col-1)*n+1):(col*n)])>seuil)
        newrez <- sample(res[((col-1)*n+1):(col*n)][abs(res[((col-1)*n+1):(col*n)])<seuil],size = length(line),replace = TRUE,prob=!is.na(res[((col-1)*n+1):(col*n)][abs(res[((col-1)*n+1):(col*n)])<seuil]))
        for(j in 1:length(line)){
          res[col,j] <- newrez[j]
        }
      }
      return(res)

    },residuals,excluded_positions,SIMPLIFY = FALSE)
  }



  samples          <- .sampling(residuals[[1]],B) # same resampling for all N triangles

  sampledResiduals <- lapply(residuals,function(x){lapply(samples,.applysample,x)})
  DFIndivBoot <- mapply(function(rez,DF,sigma,triangle){
    lapply(rez,function(.x){
      t(t(.x*sqrt(t(c(sigma^2,0)/t(triangle))))+DF)
    })},sampledResiduals,DF,sigma,triangles,SIMPLIFY = FALSE) # individual developpement factors.

  DFBoot <- mapply(function(triangle,DFind){lapply(DFind,.DFPond,triangle)},triangles,DFIndivBoot,SIMPLIFY = FALSE) # DFs



  # simulation of correlated new residuals and formating of thoose residuals.
  if(distNy == "normal"){
    # simu is a list of (n-1) B*N matrix
    simu <- lapply(1:(n-1),function(i){rmvnorm(B,sigma=rho[i,,])})
    simulation <- DFBoot
    for(c in 1:N){for(b in 1:B){
      simulation[[c]][[b]] <- sapply(simu,function(.x){.x[b,c]})
    }}

  }
  if(distNy %in% c("residuals.global","residuals")){
    samples <- .sampling(residuals[[1]],B)
    samples <- lapply(samples,function(x){as.vector(as.matrix(x)[!is.na(as.matrix(x))])[1:(n-1)]})
    simulation <- lapply(residuals,function(r){lapply(samples,function(s){as.vector(as.matrix(r))[s]})})
  }
  if(distNy == "residuals.bycolumn"){

    samples <- lapply(1:(n-1),function(i){sample(1:(n-i),size = B,replace=TRUE)})
    simulation <- DFBoot
    for(a in 1:N){for(b in 1:B){
      simulation[[a]][[b]] <- sapply(1:(n-1),function(c){residuals[[a]][c,samples[[c]][b]]})
    }}
  }




  # Calculation of next year payments on thoose residuals.
  NyCum <- mapply(function(DFBoot,diag,sigma,simulation){
    sd = sqrt((diag * c(sigma^2,0))[1:(n-1)])
    mapply(function(.x,.y){
      unname(.y * sd +  (diag * .x)[1:(n-1)])
    },DFBoot,simulation,SIMPLIFY = FALSE)},DFBoot,diag,sigma,simulation,SIMPLIFY = FALSE)


  # Coresponding incremental paymnts
  NyInc <- mapply(function(NyCum,diag){
    lapply(NyCum,function(.x){
      .x - diag[1:(n-1)]
    })},NyCum,diag,SIMPLIFY = FALSE)

  # Coreponding individual DF's
  NyDFIndiv <- mapply(function(NyCum,diag){
                        lapply(NyCum,function(.x){.x/diag[1:(n-1)]})
                      },NyCum,diag,SIMPLIFY = FALSE)


  # Coresponding Df's
  NyDF <- mapply(function(NyDFIndiv,DFIndiv,triangle){
    lapply(NyDFIndiv,.newDFPond,DFIndiv,triangle)
    },NyDFIndiv,DFIndiv,triangles,SIMPLIFY = FALSE)

  # Coresponding Ultimates
  NyUltimates <- mapply(function(triangle,NyCum,NyDF){
    mapply(function(.x,.y){
      c(.x*rev(cumprod(rev(.y[2:n]))),triangle[1,n])
    },NyCum,NyDF,SIMPLIFY = FALSE)},triangles,NyCum,NyDF,SIMPLIFY = FALSE)

  # coreponding IBNRs :
  NyIBNR <- mapply(function(y,z){
    lapply(y,function(.x){(.x-z)})
  },NyUltimates,diag,SIMPLIFY = FALSE)


  # Formating of outputs :
  # names of triangles :
  names = names(triangles)
  if(is.null(names)){
    names <- paste0("tri",1:N)
  }
  # the output will contains the marginals BootChainLadder models, retrieving the corelations will be done in external functions.
  rez <- mapply(.formatOutput,n,triangles,diag,DF,sigma,residuals,DFBoot,Ultimates,ibnr,NyCum,NyInc,NyDF,NyUltimates,NyIBNR,SIMPLIFY = FALSE)
  names(rez) <- names
  class(rez) <- c("MultiBootMackChainLadder",class(rez))
  return(rez)

}
#' mean.MultiBootMackChainLadder
#'
#' @param x A MultiBootMackChainLadder Object
#'
#' @return Data frames containing mean informations.
#' @export
#'
#' @seealso MultiBootMackChainLadder, mean.BootMackChainLadder
#'
#' @examples
#' data(ABC)
#' triangles <- list(tri1 = ABC, tri2 = ABC, tri3 = ABC)
#' MBMCL <- MultiBootmackChainLadder(triangles,100)
#' mean(MBMCL)
mean.MultiBootMackChainLadder <- function(x,...){
    margins <- lapply(x,mean)

    ByOrigin <- lapply(margins,`[[`,"ByOrigin")
    Latest <- lapply(ByOrigin,`[[`,"Latest")
    Ultimates <- lapply(ByOrigin,`[[`,"Ultimates")

    names(Latest) <- paste0("Latest.",names(Latest))
    names(Ultimates) <- paste0("Ultimates.",names(Ultimates))

    Latest <- as.data.frame(Latest)
    Ultimates <- as.data.frame(Ultimates)

    Latest$Latest.Tot = rowSums(Latest)
    Ultimates$Ultimates.Tot = rowSums(Ultimates)

    ByOr <- cbind(Latest,Ultimates)
    row.names(ByOr) <- row.names(ByOrigin[[1]])

    Totals <- as.data.frame(lapply(margins,`[[`,"Totals"))
    Totals$Tot <- rowSums(Totals)
    return(list(ByOrigin = ByOr,Totals = Totals))
}
#' CDR.MultiBootMackChainLadder
#'
#' @param MBMCL A MultiBootMackChainLadder Object
#'
#' @return Informations about CDR and IBNRS.
#' @export
#'
#' @seealso MultiBootMackChainLadder, CDR.BootMackChainLadder
#'
#' @examples
#' data(ABC)
#' triangles <- list(tri1 = ABC, tri2 = ABC, tri3 = ABC)
#' MBMCL <- MultiBootmackChainLadder(triangles,100)
#' CDR(MBMCL)
CDR.MultiBootMackChainLadder     <- function(MBMCL){

  IBNR <- as.data.frame(lapply(MBMCL,`[[`,"IBNR"))
  names(IBNR) <- paste0("IBNR.",names(IBNR))
  IBNR$IBNR.Tot = rowSums(IBNR)


  NyIBNR <- lapply(MBMCL,`[[`,"NyIBNR")
  NyIBNR$Tot <- Reduce("+", NyIBNR)
  CDR <- as.data.frame(lapply(NyIBNR,function(x){apply(x,2,sd)}))
  names(CDR) <- paste0("CDR(1).SE.",names(CDR))

  IBNR.tot <- as.vector(colSums(IBNR))
  CDR.tot = lapply(NyIBNR,rowSums)
  CDR.tot = as.data.frame(lapply(CDR.tot,sd))


  ByOrigin = cbind(IBNR,CDR)
  Totals = data.frame(IBNR = IBNR.tot,`CDR(1)S.E.` = t(CDR.tot))

  return(list(ByOrigin = ByOrigin,Totals=Totals))
}
#' Title
#'
#' @param MBMCL A MultiBootMackChainLadder Object
#' @param ByOrigin If Set to TRUE, Next year IBNRS will be rgiveng with per-origin details.
#'
#' @return Data frame with next year informations.
#' @export
#'
#' @examples
#' data(ABC)
#' triangles <- list(tri1 = ABC, tri2 = ABC, tri3 = ABC)
#' MBMCL <- MultiBootmackChainLadder(triangles,100)
#' NyIBNR(MBMCL)
NyIBNR <- function(MBMCL,ByOrigin = FALSE){
  NyIBNR <- lapply(MBMCL,`[[`,"NyIBNR")
  NyIBNR$Tot <- Reduce("+", NyIBNR)
  if(!ByOrigin){
    NyIBNR = lapply(NyIBNR,rowSums)
  }
  return(as.data.frame(NyIBNR))
}
#' Title
#'
#' @param MBMCL A MultiBootMackChainLadder Object
#' @param ... Elements to be passed to the cor function
#'
#' @return A corelation matrix of next-year IBNRs.
#' @export
#'
#' @examples
#' data(ABC)
#' triangles <- list(tri1 = ABC, tri2 = ABC, tri3 = ABC)
#' MBMCL <- MultiBootmackChainLadder(triangles,100)
#' Corel(MBMCL)
Corel <- function(MBMCL,...){
  NyIBNR <- NyIBNR(MBMCL)
  NyIBNR$Tot <- NULL
  return(cor(NyIBNR,...))
}


##### Specific structure for my memoire                                              #####

#' Analyse
#'
#' This function perform the analysis from my memoire on a set of triangles.
#'
#'
#' @param triangles the set of triangles. Should be a named list with 3 triangles for PSAPs and the last one for PSNEMs.
#' @param B Number of bootstrap replicates
#' @param distNy1 Parameter for the MultiBookMackChainLadder part (psaps) of the procedure.
#' @param distNy2 Parameter for the BookMackChainLadder part (psnems) of the procedure.
#' @param seuil Parameter for both parts.
#' @param pdd Should be a data.frame with one row and 3 columns.
#' @param Dossiers Working folders informaitons.
#' @param Zonnage.capi Should the capitalisation bootstrap be zonned ?
#'
#' @return an "Etudeprincipale" object wich has only a print method and contains everything (try to str it !)
#' @export
Analyse <- function(triangles,B=100,distNy1="normal",distNy2="normal",seuil=NA,pdd = NA, Dossiers = NA,Zonnage.capi=FALSE){

  # Achtung : Le 4eme triangle doit etre celui en capi.
  # Dans les 3 premiers triangles ( ils doivent etre nomes), l'un d'entre eux doit s'apeller RCDO.
  cat("\n MBMCL (1/4)...\n")
  # Calcul du bootstrapp syncrho sur les triangles SI :
  MBMCL <- MultiBootMackChainLadder(triangles[1:3],
                                    B=B,
                                    names=names(triangles)[1:3],
                                    distNy = distNy1,
                                    seuil=seuil)

  cat("tri.ds.ult (2/4)... \n")
  # Passage a l'utlime du triangle de RCDO en droc-surv sur chaque nouvelle cadence :
  tri.ds.ult <- MBMCL$RCDO$NyUltimates %>%
    apply(1,list) %>%
    unlist(recursive=FALSE) %>%
    map2(list(.diag(triangles$RCDO)),~.y/.x) %>%
    map(~incr2cum(.passage.a.lultime(triangles[[4]],.x)))

  cat("BMCL.capis (3/4)... \n")
  # Bootstrapp du triangle RCDO ultime :
  BMCL.capis <- tri.ds.ult %>% map(BootMackChainLadder,B=1,distNy = distNy2,seuil=seuil,zonnage=Zonnage.capi)

  class(BMCL.capis) <- c("BootMackChainLadderCapi",class(BMCL.capis))


  valeur.psnem.deterministe <- (.diag(triangles$RCDO)/.ultimates(triangles$RCDO)) %>%
  {.passage.a.lultime(triangles[[4]],.)} %>%
    incr2cum %>%
    .ibnr %>%
    sum
  valeur.psap.capi.deterministe <- (.diag(triangles$RCDO)/.ultimates(triangles$RCDO)) %>%
  {sum(.diag(incr2cum(.passage.a.lultime(triangles[[4]],.)))) - sum(triangles[[4]],na.rm=TRUE)}


  # recuperation des 5 provisions uniquement pour analyse :
  ibnr <- MBMCL %>%
    map("NyIBNR") %>%
    map(rowSums) %>%
    enframe %>%
    add_row(name = "RCDO.psnem",value = list(BMCL.capis %>% map("NyIBNR") %>% map_dbl(sum))) %>%
    spread(name,value) %>%
    unnest %>%
    mutate(RCDO.ibnr.capi = (BMCL.capis %>% map("Latest") %>% map_dbl(sum)) - sum(triangles[[4]],na.rm=TRUE))

  psap <- ibnr %>% transmute(
    RCDC = RCDC + pdd$RCDC,
    RCDO = RCDO.ibnr.capi + pdd$RCDO, ######### On pourrais prendre ici les psap calcul?es en repart, mais on prend celles en capi.
    RCG = RCG + pdd$RCG,
    PSNEM = RCDO.psnem,
    LOB8 = RCDC + RCDO + PSNEM + RCG
  )

  # petite mise en jambe sur les residus :
  analyse.residus <- map(triangles[1:3],.residuals) %>%
    map(as.data.frame) %>%
    map(as.tibble) %>%
    map2(names(.),~mutate(.x,tri = .y)) %>%
    {do.call(rbind,.)} %>%
    select(tri,origin,dev,value) %>%
    mutate(is.ok = map_dbl(value, function(x){return(x>=-2 & x <= 2)}))

  cat("out (4/4)... \n")
  out <- list(triangles = triangles,
              B = B,
              distNy1 = distNy1,
              distNy2 = distNy2,
              seuil = seuil,
              pdd = pdd,
              Dossiers = Dossiers,
              MBMCL = MBMCL,
              BMCL.capis = BMCL.capis,
              valeur.psnem.deterministe = valeur.psnem.deterministe,
              valeur.psap.capi.deterministe = valeur.psap.capi.deterministe,
              ibnr = ibnr,
              psap = psap,
              analyse.residus = analyse.residus
              )
  class(out) <- c("Etudeprincipale",class(out))
  return(out)

}

#' print.Etudeprincipale
#'
#' @param x An EtudePrincipale Object
#'
#' @return returns x after aving printed a lot of informations about it.
#' @export
#'
print.Etudeprincipale <- function(x,...){
  cat("Ceci est une etude bootstrap synchro en LOB8. \n\n1. Parametres :\n \tB =",x$B,
      "\n\tdistNy1 = ",x$distNy1,
      "\n\tdistNy2 = ",x$distNy2)
  if(!is.na(x$seuil)){
  cat("\n\tSeuil limite des residus = ",x$seuil)
  }
  cat("\n\tPDD Totale :  = ",sum(x$pdd))
  cat("\n\t Taille des triangles : ",dim(x$triangles[[1]])[[1]])
  # cat("\n\nDossiers d'enregistrement: ")
  # print(rez$Dossiers %>% unlist)

  # petit check : les ultimes de l'annee prochaine moins les ibnr de l'annee prochaine moins la diag de cette annee, ca doit faire 0.
  if(c((rowSums(x$MBMCL$RCDO$NyUltimates) - rowSums(x$MBMCL$RCDO$NyIBNR) - sum(.diag(x$triangles$RCDO))),
    (rowSums(x$MBMCL$RCDC$NyUltimates) - rowSums(x$MBMCL$RCDC$NyIBNR) - sum(.diag(x$triangles$RCDC))),
    (rowSums(x$MBMCL$RCG$NyUltimates)  - rowSums(x$MBMCL$RCG$NyIBNR)  - sum(.diag(x$triangles$RCG)))) %>% round(6) %>% unique != 0){
    cat("Un des tests echous, le modele possede probablement une erreur de precision quelque part.")
  }

  cat("\n\n2. Totaux du bootstrap synchro partie repartition :\n")
  print(mean(x$MBMCL)$Totals)

  cat("\n\n3. PSNEM et PSAP deterministes obtenus en capi :\n")
  cat("\n\tPSNEM = ",x$valeur.psnem.deterministe)
  cat("\n\tPSAP  = ",x$valeur.psap.capi.deterministe)


  cat("\n\n4. IBNR moyens (EN RCDO : calcule en repart ou en capi) :\n")
  print(colMeans(x$ibnr))
  cat("\nCorrelation des IBNR : \n")
  print(cor(x$ibnr))

  cat("\n\n5. PSAP (en moyenne) :\n")
  print(colMeans(x$psap))
  cat("\nCorrelation des PSAP : \n")
  print(cor(x$ibnr))


  cat("\n\n6. Resultats en termes de sigma de reserve :\n")
  x$psap %>%
    map_dfr(~data.frame(sd = sd(.),mean = mean(.))) %>%
    mutate(grh = colnames(x$psap)) %>%
    select(grh,sd,mean) %>%
    mutate(sigma = sd/mean) %T>% print


  cat("\n\n7. Pourcentage de residus en dehors de (-2;2) : ",table(x$analyse.residus$is.ok) %>% {.[1]/(.[1]+.[2])})
  return(NULL)
}

#' plots of EtudePrincipale
#'
#'  A serie of plotting functions exists for an etudeprincipale.
#'  Some more will probably be added later.
#'
#'
#' @param x A EtudePrincipale object
#'
#' @return dependns. sometimes thees functions returns a plot.
#' @export
#'
#' @details Following plotting functions are avaliable : plot_factors_density plot_reserve_density plot_ibnr_density plot_burn_in plot_rho plot_resid_norm plot_resid_dens plot_resid_dens2 plot_resid_stabilite plot_resid_margins plot_cdr_mw plot_cadences_capi
#'
#' @aliases plot_factors_density plot_reserve_density plot_ibnr_density plot_burn_in plot_rho plot_resid_norm plot_resid_dens plot_resid_dens2 plot_resid_stabilite plot_resid_margins plot_cdr_mw plot_cadences_capi
plot_factors_density <- function(x){

  # graphique des cadences de developpement bootstrapees
    x$MBMCL %>%
    map("DFBoot") %>%
    map(as.data.frame) %>%
    map(gather) %>%
    enframe %>%
    unnest %>%
    set_names(c("grh","dev","value")) %>%
    mutate(dev = as.numeric(dev)) %>%
    group_by(grh,dev) %>%
    summarise(mean = mean(value),
              lower=quantile(value,0.01),
              upper=quantile(value,0.99)) %>% ungroup %>%
    ggplot(aes(x = dev,y = mean,color=grh)) +
    facet_grid(grh ~ .) +
    geom_line(size=1,color="black") +
    geom_ribbon(aes(ymin=lower,ymax=upper,fill=grh),size=0, alpha=0.5)+
    geom_hline(yintercept = 1,size=0.5) +
    facet_grid(grh ~ .,scales = "free") +
    theme(legend.position="none") +
    ggtitle("Facteurs de developpements bootstrappes en repartition")

}
#' @export
plot_reserve_density <- function(x){
  p2 <- x$psap %>%
    gather %>%
    ggplot(aes(x=value,col=key,fill=key)) +
    geom_density(alpha=0.1) +
    ggtitle("Densitees des provisions des differents triangles")

  # densite des psap et des psnem :
  p3 <- data.frame(psnem = x$psap$PSNEM, psap = x$psap$RCDO) %>%
    gather %>%
    ggplot(aes(x=value,col=key,fill=key)) +
    geom_density(alpha=0.1) +
    ggtitle("Densitees des psap et des psnem (a un an) sur la RCDO")


  # densiite des IBNR :
  p4 <- x$ibnr %>%
    gather %>%
    ggplot(aes(x=value,col=key,fill=key)) +
    geom_density(alpha=0.1) +
    ggtitle("Densitees des ibnr uniquement")


  grid.arrange(grobs=list(p2,p3,p4),layout_matrix = rbind(c(1),c(2),c(3)))
  # rm(p1,p2,p3,p4)
  # plot(1)
  # plot(p5)
}

#' @export
plot_ibnr_density<- function(x){

  par(mfrow=c(3,2))

  x$ibnr$RCG %>% density %>% plot(main="IBNR en RCG")
  abline(v = x$triangles[["RCG"]] %>% .ibnr %>% sum)

  x$ibnr$RCDC %>% density %>% plot(main="IBNR en RCDC")
  abline(v = x$triangles[["RCDC"]] %>% .ibnr %>% sum)

  x$ibnr$RCDO %>% density %>% plot(main="IBNR en RCDO (repart)")
  abline(v = x$triangles[["RCDO"]] %>% .ibnr %>% sum)

  x$ibnr$RCDO.ibnr.capi %>% density %>% plot(main="IBNR en RCDO (capi)")
  abline(v = x$valeur.psap.capi.deterministe)

  x$ibnr$RCDO.psnem %>% density %>% plot(main="PSNEM RCDO")
  abline(v=x$valeur.psnem.deterministe) # on a donc bien une psnem centre sur la psnem deterministe.

  plot(density((x$ibnr$RCDO.ibnr.capi-x$ibnr$RCDO)/x$ibnr$RCDO), main="Erreur relative des ibnr RCDO capi/repart")
  abline(v = mean((x$ibnr$RCDO.ibnr.capi-x$ibnr$RCDO)/x$ibnr$RCDO))
  par(mfrow = c(1,1))
}
#' @export
plot_burn_in <- function(x){
  # Burn-in de la varianc des provisions :
  p1 <- x$psap %>%
    mutate_all(.cumsd) %>%
    mutate(row = row_number()) %>%
    gather(key="key",value="value",-row) %>%
    ggplot(aes(x= row, y=value,col=key)) +
    geom_line() +
    ggtitle("Burn-in de la variance des provisions")

  p2 <- x$psap %>%
    mutate_all(~.cumsd(.x)/cummean(.x)) %>%
    mutate(row = row_number()) %>%
    gather(key="key",value="value",-row) %>%
    ggplot(aes(x= row, y=value,col=key)) +
    geom_line() +
    ggtitle("Burn-in du sigma de reserve")

  p3 <- x$psap %>% select(LOB8) %>%
    mutate_all(~.cumsd(.x)/cummean(.x)) %>%
    mutate(row = row_number()) %>%
    gather(key="key",value="value",-row) %>%
    ggplot(aes(x= row, y=value,col=key)) +
    geom_line() +
    ggtitle("Burn-in du sigma de reserve LOB")

  p4 <- x$psap %>%
    .cumcorel %>% select(2:5,8:10,14:15,20) %>%
    mutate(row = row_number()) %>%
    gather(key="key",value="value",-row) %>%
    ggplot(aes(x= row, y=value,col=key)) +
    geom_line() +
    ggtitle("Burn-in des corelations")

  grid.arrange(p1,p2,p3,p4)
}
#' @export
plot_rho <- function(x){
  # graphique des coeficients de corelations des triangles.
  rho = .rho(x$triangles)
  map(1:(as.numeric(rev(rownames(x$triangles[[1]]))[1])-1990-1),~c(rho[.x,,][1,2:3],rho[.x,,][2,3])) %>% map_dbl(1) %>% plot(col=4,main="Parametres de corelations des triangles.",ylim=c(-1,1))
  map(1:(as.numeric(rev(rownames(x$triangles[[1]]))[1])-1990-1),~c(rho[.x,,][1,2:3],rho[.x,,][2,3])) %>% map_dbl(2) %>% points(col=2)
  map(1:(as.numeric(rev(rownames(x$triangles[[1]]))[1])-1990-1),~c(rho[.x,,][1,2:3],rho[.x,,][2,3])) %>% map_dbl(3) %>% points(col=3)
  abline(h=0)
}
#' @export
plot_resid_norm <- function(x){
  par(mfrow = c(2,3))
  x$analyse.residus %>%
    group_by(tri) %>%
    summarise(residus = list(value)) %>%
    {map2(.$residus,.$tri,function(rez,nom){
      qqnorm(rez,main=paste0("QQ-Norm : ",nom))
      abline(0,1)
    })}

  x$analyse.residus %>% filter(tri == "RCDO") %>% .$value %>% density(na.rm=TRUE) %>% plot(xlim=c(-5,5),main="Densites des residus en RCDO")
  curve(dnorm,add=TRUE,col=2)
  x$analyse.residus %>% filter(tri == "RCDC") %>% .$value %>% density(na.rm=TRUE) %>% plot(xlim=c(-5,5),main="Densites des residus en RCDC")
  curve(dnorm,add=TRUE,col=2)
  x$analyse.residus %>% filter(tri == "RCG") %>% .$value %>% density(na.rm=TRUE) %>% plot(xlim=c(-5,5),main="Densites des residus en RCG")
  curve(dnorm,add=TRUE,col=2)
  par(mfrow = c(1,1))
}
#' @export
plot_resid_dens <- function(x){
  # graphiques des residus
  x$analyse.residus %>%
    drop_na %>%
    select(tri,value) %>%
    ggplot(aes(x = value,color=tri,fill=tri)) +
    geom_density(alpha = 0.1) +
    ggtitle("Densite des residus des triangles en survenance-inventaire")
}
#' @export
plot_resid_dens2 <- function(x){
  # analyse des residus :

  rez <- x$analyse.residus %>% drop_na %>% .$value
  hist(rez,freq=FALSE,main="Distribution des residus / densitee d'une N(0,1)")
  x = seq(min(rez),max(rez),length=1000)
  y = dnorm(x)
  lines(x,y)
}
#' @export
plot_resid_stabilite <- function(x){

  p1 <- x$analyse.residus %>% drop_na %>%
    mutate(dev = as.numeric(dev)) %>%
    ggplot(aes(x=dev,y=value,color=tri)) +
    geom_point() +
    geom_smooth(method="loess") +
    ggtitle("Residus des differents triangles par annee de developpement")
  p2 <- x$analyse.residus %>% drop_na %>%
    mutate(origin = as.numeric(origin)) %>%
    ggplot(aes(x=origin,y=value,color=tri)) +
    geom_point() +
    geom_smooth(method="loess") +
    ggtitle("Residus des differents triangles par annee d'origine")
  p3 <- x$analyse.residus %>%
    mutate(cal = as.numeric(origin)+as.numeric(dev)) %>% drop_na %>%
    ggplot(aes(x=cal,y=value,color=tri)) +
    geom_point() +
    geom_smooth(method="loess") +
    ggtitle("Residus des differents triangles par annee calendaire")

  grid.arrange(p1,p2,p3)
}
#' @export
plot_resid_margins <- function(x){
  # et finalement les planches de graphiques de mack-chain-ladder
  suppressWarnings(x$triangles[1:3] %>% map(MackChainLadder,est.sigma = "Mack") %>% map(plot))
  suppressWarnings(x$triangles[[4]] %>% incr2cum %>% MackChainLadder(est.sigma = "Mack") %>% plot)
}
#' @export
plot_cdr_mw <- function(x){
  # petite comparaison des CDR avec ceux implementes dans la librairie ChainLadder :
  cdr.booted <- c(CDR(x$MBMCL$RCDO)$`CDR(1)S.E.`[1:24],CDR(x$MBMCL$RCDC)$`CDR(1)S.E.`[1:24],CDR(x$MBMCL$RCG)$`CDR(1)S.E.`[1:24])
  suppressWarnings(classical.cdr.merz <- c(x$triangles$RCDO %>% MackChainLadder %>% CDR %>% .$`CDR(1)S.E.` %>% .[1:24],
                                           x$triangles$RCDC %>% MackChainLadder %>% CDR %>% .$`CDR(1)S.E.` %>% .[1:24],
                                           x$triangles$RCG %>% MackChainLadder %>% CDR %>% .$`CDR(1)S.E.` %>% .[1:24]))
  plot(cdr.booted,classical.cdr.merz,main="CDR(1) : bootstrap VS M&W classique")
  abline(0,1)
}
#' @export
plot_cadences_capi <- function(x,type =2){


  DFtoCad <- . %>%
    c(1,.) %>%
    cumprod %>%
    {.-c(1,.)} %>%
    {c(1,.[2:(length(.)-1)])} %>%
    {./sum(.)} %>% set_names(paste0("cad",1:length(.)))


  cadences_boot <- x$BMCL.capis %>% map("DF") %>% map(DFtoCad)

  if(type==1){
   p <-  cadences_boot %>%
      unlist %>%
      {list(.,names(.))} %>%
      enframe %>%
      spread(key = "name",value = "value") %>%
      unnest %>%
      set_names(c("cad","date")) %>%
      select(date,cad) %>%
      separate(date,c("plop","date"),sep = "d") %>%
      transmute(x = as.numeric(date),y=cad) %>%
      #filter(x < 19) %>%
      mutate(
        zonne = as.factor(1*(x<3)+2*(x>2)*(x<10)+3*(x>9)*(x<13)+4*(x>12))
      ) %>%
      ggplot(aes(as.factor(x),y,color=zonne))+
      ggtitle("Cadence de manifestations de la RCDO - Zonnee.") +
      geom_boxplot(width=1,outlier.size = 0.7)+
     stat_summary(fun.y=mean, geom="line", aes(group=1))
    #stat_summary(fun.y=mean, geom="line", aes(group=1))  +
    #stat_summary(fun.y=mean, geom="point")# +
    #geom_vline(xintercept = c(2.5,9.5,12.5))
  } else {
    p <- cadences_boot %>%
      unlist %>%
      {list(.,names(.))} %>%
      enframe %>%
      spread(key = "name",value = "value") %>%
      unnest %>%
      set_names(c("cad","date")) %>%
      select(date,cad) %>%
      separate(date,c("plop","date"),sep = "d") %>%
      transmute(x = as.numeric(date),y=cad) %>%
      filter(x < 19) %>%
      mutate(
        zonne = as.factor(1*(x<3)+2*(x>2)*(x<10)+3*(x>9)*(x<13)+4*(x>12))
      ) %>%
      ggplot(aes(as.factor(x),y,color=zonne))+
      ggtitle("Cadence de manifestations de la RCDO - Zonnee.") +
      geom_boxplot(width=1,outlier.size = 0.7)+
      #stat_summary(fun.y=mean, geom="line", aes(group=1))  +
      #stat_summary(fun.y=mean, geom="point") +
      geom_vline(xintercept = c(2.5,9.5,12.5)) +
      stat_summary(fun.y=mean, geom="line", aes(group=1))
  }

  return(p)
}



#' .creer_nom_fichier
#'
#' @param graves TRue or false
#' @param seuil numeric or NA
#' @param type residuals or normal
#' @param data Reglements or Charges
#' @param annee 2014, 2015, 2016 and 2017 are availiable
#' @param zonnage True or False (zonning the capitalisation residuals)
#'
#' @return The file name (a string)
#' @export
.creer_nom_fichier <- function(graves,seuil,type,data,annee,zonnage){
  # enplacement du fichier de sauvegarde :
  if(graves){
    .grv = "avecGraves"
  } else {
    .grv = "sansGraves"
  }
  if(!is.na(seuil)){
    .seuil = paste0("SeuilDe",seuil)
  } else {
    .seuil = "sansSeuil"
  }
  if(type == "residuals"){
    .resid = "typeRezi"
  } else {
    if(type == "normal"){
      .resid = "typeNorn"
    }
  }
  if(data == "Reglements"){
    .data = "Reg"
  } else {
    if(data == "Charges"){
      .data = "Chg"
    }
  }
  if(zonnage){
    .zon = "avecZonne"
  } else {
    .zon = "sansZonne"
  }
  .fichier.resultats = paste0(paste(.data,annee,.resid,.grv,.seuil,.zon,sep="_"),".Rdata")
  return(.fichier.resultats)
}

##### Tests                                                                          #####
#
#
# # test :
# library(ChainLadder)
# library(magrittr)
#
#
#
# data(ABC)
# data(auto)
# MCL <- MackChainLadder(ABC,est.sigma="Mack")
# BMCL1 <- BootMackChainLadder(ABC,B=10000,distNy = "normal")
# BMCL2 <- BootMackChainLadder(ABC,B=10000,distNy = "residuals")
#
# # r?sum? :
# BMCL1
# BMCL2
# MCL
#
# # cdr :
# plot(as.numeric(rownames(CDR(MCL))[1:11]),CDR(MCL)$`CDR(1)S.E.`[1:11],type="l")
# points(as.numeric(rownames(CDR(MCL))[1:11]),CDR(BMCL1)$`CDR(1)S.E.`[1:11],type="l",col=2)
# points(as.numeric(rownames(CDR(MCL))[1:11]),CDR(BMCL2)$`CDR(1)S.E.`[1:11],type="l",col=3)
#
# # the convergence is neat :)
#
#
# # lets try on other triangles :
# data(MW2008)
#
# MCL <- MackChainLadder(MW2008,est.sigma="Mack")
# BMCL1 <- BootMackChainLadder(MW2008,B=10000,distNy = "normal")
# BMCL2 <- BootMackChainLadder(MW2008,B=10000,distNy = "residuals")
#
# # r?sum? :
# BMCL1
# BMCL2
# MCL
#
# # cdr :
# plot(as.numeric(rownames(CDR(MCL))[1:11]),CDR(MCL)$`CDR(1)S.E.`[1:11],type="l")
# points(as.numeric(rownames(CDR(MCL))[1:11]),CDR(BMCL1)$`CDR(1)S.E.`[1:11],type="l",col=2)
# points(as.numeric(rownames(CDR(MCL))[1:11]),CDR(BMCL2)$`CDR(1)S.E.`[1:11],type="l",col=3)
#
#
# MBMCL1 <- MultiBootMackChainLadder(auto,B=10000,distNy = "normal")
# MBMCL2 <- MultiBootMackChainLadder(auto,B=10000,distNy = "residuals.global")
# MBMCL3 <- MultiBootMackChainLadder(auto,B=10000,distNy = "residuals.bycolumn")
# MBMCL1
# MBMCL2
# MBMCL3
# data.frame(cdr1 = CDR(MBMCL1)$ByOrigin$`CDR(1).SE.Tot`,
#            cdr2 = CDR(MBMCL2)$ByOrigin$`CDR(1).SE.Tot`,
#            cdr3 = CDR(MBMCL3)$ByOrigin$`CDR(1).SE.Tot`) %>% matplot(type="l")
# legend("topleft",cex=0.7,col=c(1,2,3),lty=c(1,2,3),legend = c("H : residus normaux, correles, simules","H : residus globalement correles, boostrapes","H : residus correles par colonne, bootstrape"))
#
#
# # sigma de reserve total :
# sigma1 <- CDR(MBMCL1)$Total["Tot",2]/CDR(MBMCL1)$Total["Tot",1]
# sigma2 <- CDR(MBMCL2)$Total["Tot",2]/CDR(MBMCL2)$Total["Tot",1]
# sigma3 <- CDR(MBMCL3)$Total["Tot",2]/CDR(MBMCL3)$Total["Tot",1]
#

##### Fin #####
