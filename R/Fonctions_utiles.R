################################################################################
################################# Fonctions utilitaires.
################################################################################
.check.triangle         <- function(triangle,dim=dim(triangle)){
  # FOnction de vÃ©rification des dimentions d'un triangle :
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
    # Passage en incrémental si et seulement si ça renvois pas une erreur :
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

################################################################################
################################# Passage a l'ultime
################################################################################
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

################################################################################
################################# Fonction de calcul cumulé
################################################################################
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


################################################################################
################################# Création des noms des fichiers
################################################################################
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
