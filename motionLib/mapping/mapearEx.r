library(sp)
library(RColorBrewer)
library(colorspace)
library(ggplot2)
library(rgeos)
library(maptools)
library(directlabels)
library(ggrepel)
library(Cairo)
library(mapproj)
paletasInvertidas <- c('Spectral', 'RdBu', 'RdYlBu', 'RdYlGn')

crearEscala <- function(escala, colores=NULL, brewerPal='Spectral', 
                        invertirPaleta=brewerPal %in% paletasInvertidas, continuo=F,
                        iniciosIntervalosIsoLineas=NULL) {
  u <- unique(escala)
  
  if (is.null(colores)) {
    # Si la escala es continua, los colores son de los puntos del intervalo,
    # si es discreta, los colores son de los intervalos, y por lo tanto hay uno menos
    if (continuo) { nColores <- length(u)
    } else { nColores <- length(u) - 1 }
    
    if (nColores < 3) {
      colores <- brewer.pal(3, name = brewerPal)
      if (nColores == 1) { colores <- colores[1]
      } else { colores <- colores[c(1,3)]}
    } else if (nColores <= brewer.pal.info[brewerPal,]$maxcolors) {
      colores <- brewer.pal(nColores, name = brewerPal)
    } else {
      f <- colorRampPalette(brewer.pal(brewer.pal.info[brewerPal,]$maxcolors, name = brewerPal))
      colores <- f(nColores) 
    }
    if (invertirPaleta) colores <- rev(colores)
  }
  
  if (length(u) < length(escala)) {
    i <-  match(u, escala)
    escala <- escala[i]
    
    if (length(escala) == 1) {
      escala <- c(escala, escala[1] + 1)
      colores <- colores[c(1, length(colores))]
    } else if (i[length(i)] > length(colores)) { 
      colores <- colores[i[1:(length(i)-1)]]
    } else { 
      colores <- colores[i] 
    }
  }
 
  if (length(colores) == length(escala) - 1 && !continuo) colores <- c(colores, '')
  return (list(escala=escala, colores=colores, continuo=continuo, iniciosIntervalosIsoLineas=iniciosIntervalosIsoLineas))
}

crearEscalaTemperaturaWRF <- function() {
  colR <- c(40,60,30,0,75,200,200,160,60,0,0,255,255,240,255,255,230,200,150,110,70,255,255,250,205,130,220,160,72)
  colG <- c(0,40,60,150,240,255,240,230,255,210,200,255,220,130,70,0,200,150,80,40,0,180,80,0,0,0,220,140,60)
  colB <- c(160,180,255,255,255,255,150,50,30,140,0,170,0,40,0,0,45,0,0,0,0,255,250,180,205,220,255,255,200)
  escala <- seq.int(from = -12, to = 46, by = 2)
  colores <- rgb(colR, colG, colB, 255, maxColorValue = 255)
  return(crearEscala(escala, colores))
}

achicarToNDigitos <- function(x, nDigitos=1) {
  pot <- 10^-nDigitos
  return(floor(x / pot) * pot)
}

agrandarToNDigitos <- function(x, nDigitos=1) {
  pot <- 10^-nDigitos
  return(ceiling(x / pot) * pot)
}

rangoExtendidoANDigitos <- function(x, nDigitos=1, na.rm=T) {
  rango <- range(x, na.rm = na.rm)
  rango[1] <- achicarToNDigitos(rango[1], nDigitos = nDigitos)
  rango[2] <- agrandarToNDigitos(rango[2], nDigitos = nDigitos)
  return(rango)
}

crearEscalaDosPuntos <- function(inicio, fin, colores=NULL, nIntervalos=9, intervaloFinalCerrado = TRUE, brewerPal='Spectral', invertirPaleta=brewerPal %in% paletasInvertidas,  continuo=F) {
  escala <- seq(from= inicio, to = fin, length.out = nIntervalos)
  if (is.null(colores)) {
    if (intervaloFinalCerrado && !continuo) { colores <- brewer.pal(nIntervalos - 1, brewerPal)
    } else { colores <- brewer.pal(nIntervalos, brewerPal) }
    if (invertirPaleta) colores <- rev(colores)
  }
  return(crearEscala(escala = escala, colores = colores, continuo = continuo, invertirPaleta = F))
}

crearEscalaTresPuntos <- function(inicio, medio, fin, colorInicio='blue', colorMedio='white', colorFin='red', nIntervalos=9,
                                  intervaloFinalCerrado = TRUE, space='Lab', continuo=F) {
  delta <- (fin - inicio) / (nIntervalos - 1)
  invDelta <- 1 / delta
  
  if (intervaloFinalCerrado && !continuo) { nColores <- nIntervalos - 1
  } else { nColores <- nIntervalos }
  
  iMedio <- trunc((medio - inicio) * invDelta) + 1
  
  coloresInicioMedio <- colorRampPalette(colors = c(colorInicio, colorMedio), space=space)(iMedio)
  coloresMedioFin <- colorRampPalette(colors = c(colorMedio, colorFin), space=space)(nColores - iMedio + 1)
  colores <- c(coloresInicioMedio, coloresMedioFin[2:length(coloresMedioFin)])
  
  escalaInicioMedio <- seq(from = inicio, to = medio, length.out = iMedio)
  escalaMedioFin <- seq(from = medio, to = fin, length.out = nIntervalos - iMedio + 1)
  escala <- c(escalaInicioMedio, escalaMedioFin[2:length(escalaMedioFin)])
  
  return(crearEscala(escala = escala, colores = colores, continuo = continuo))
}

crearEscalaEquiespaciada <- function(datos, nDigitos=1, nIntervalos=9, colores=NULL, brewerPal='Spectral', invertirPaleta=brewerPal %in% paletasInvertidas, continuo=F) {
  rango <- range(datos, na.rm=T)
  escala <- unique(round(seq(from=achicarToNDigitos(rango[1], nDigitos = nDigitos), to=agrandarToNDigitos(rango[2], nDigitos = nDigitos), length.out=nIntervalos+1), nDigitos))
  
  return (crearEscala(escala=escala, colores=colores, brewerPal = brewerPal, invertirPaleta=invertirPaleta, continuo=continuo))
}

crearEscalaEquiespaciadaRecortada <- function(datos, nDigitos=1, nIntervalos=9, colores=NULL, brewerPal='Spectral', invertirPaleta=brewerPal %in% paletasInvertidas,
                                              continuo=F, lowerTrimPU=0.005, upperTrimPU=lowerTrimPU) {
  datos <- sort(na.omit(datos))
  diffDatos <- diff(datos)
  
  aux <- which(diffDatos - mean(diffDatos) > sd(diffDatos))
  aux2 <- diff(aux)
  i <- which.max(aux2)
  iLowerTrim <- aux[i]
  iUpperTrim <- aux[i+1]
    
  # plot(1:(length(datos)-2), diff(diffDatos))
  # plot(1:(length(datos)-1), diffDatos)
  # plot(1:(length(datos)), datos)
  
  #modelo <- lm(formula = y~I(x^3)+I(x^2)+x+1, data = data.frame(y=datos, x=1:length(datos)))
  #plot(1:(length(datos)), datos, pch='.')
  #lines(1:(length(datos)), modelo$fitted.values, col='red')
  
  #iLowerTrim <- which(modelo$fitted.values < datos)[1]
  #aux <- which(modelo$fitted.values > datos)
  #iUpperTrim <- aux[length(aux)]
  
  #iLowerTrim <- round(length(datos) * lowerTrimPU)
  #iUpperTrim <- round(length(datos) * (1 - upperTrimPU))
  iLowerTrim <- max(iLowerTrim, 1)
  iUpperTrim <- min(iUpperTrim, length(datos))
  
  nIntervalosAux <- nIntervalos
  if (iLowerTrim > 1) nIntervalosAux <- nIntervalosAux - 1
  if (iUpperTrim < length(datos)) nIntervalosAux <- nIntervalosAux - 1
  
  rango <- c(datos[iLowerTrim], datos[iUpperTrim])
  escala <- unique(round(seq(from=achicarToNDigitos(rango[1], nDigitos = nDigitos), to=agrandarToNDigitos(rango[2], nDigitos = nDigitos), length.out=nIntervalosAux+1), nDigitos))
  
  if (iLowerTrim > 1) {
    if (iUpperTrim < length(datos)) { escala <- c(datos[1], escala, datos[length(datos)])
    } else { escala <- c(datos[1], escala) }
  } else if (iUpperTrim < length(datos)) { escala <- c(escala, datos[length(datos)])  }
  
  escala <- crearEscala(escala=escala, colores=colores, brewerPal = brewerPal, invertirPaleta=invertirPaleta, continuo=continuo)
  return (ajustarExtremosEscala(escala=escala, datos = datos, nDigitos = nDigitos))
}

crearEscalaEquiespaciadaDistinguirMenoresOIgualesAUmbral <- function(datos, umbral, nDigitos=1, nIntervalos=9, colores=NULL, brewerPal='Spectral', 
                                                                     invertirPaleta=brewerPal %in% paletasInvertidas, continuo=F) {
  # Crea una escala que tiene en el primer intervalo todos los valores menores o iguales a
  # umbral, y en el resto de los intervalos es equiespaciada desde el mínimo valor en datos mayor a umbral 
  # hasta el máximo valor
  if (any(datos > umbral, na.rm = T)) { rango <- range(datos[datos > umbral], na.rm=T)
  } else { rango <- range(datos, na.rm = T) }
  escala <- unique(c(min(datos, na.rm=T), round(seq(from=achicarToNDigitos(rango[1], nDigitos = nDigitos), to=agrandarToNDigitos(rango[2], nDigitos = nDigitos), length.out=nIntervalos-1), nDigitos)))
  
  return (crearEscala(escala=escala, colores=colores, brewerPal = brewerPal, invertirPaleta=invertirPaleta, continuo=continuo))
}

crearEscalaEquiespaciadaLogaritmica <- function(datos, nDigitos=1, nIntervalos=9, colores=NULL, brewerPal='Spectral', invertirPaleta=brewerPal %in% paletasInvertidas, continuo=F) {
  rango <- range(datos, na.rm=T)
  if (rango[1] <= 0) stop('mapearEx.crearEscalaEquiespaciadaLogaritmica: para una escala logaritmica todos los valores deben ser > 0')
  
  rango <- log(rango)
  escala <- round(seq(from=achicarToNDigitos(rango[1], nDigitos = nDigitos), to=agrandarToNDigitos(rango[2], nDigitos = nDigitos), length.out=nIntervalos+1), nDigitos)
  escala <- unique(round(exp(escala), nDigitos))
  
  return (crearEscala(escala=escala, colores=colores, brewerPal=brewerPal, continuo=continuo))
}

crearEscalaEquiespaciadaLogaritmicaDistinguirMenoresAUmbral <- function(datos, umbral, base=exp(1), nDigitos=1, nIntervalos=9, colores=NULL, brewerPal='Spectral', invertirPaleta=brewerPal %in% paletasInvertidas, continuo=F) {
  #datos <- observaciones$R3
  if (umbral <= 0) stop('mapearEx.crearEscalaEquiespaciadaLogaritmica: para una escala logaritmica todos los valores deben ser > 0')
  
  rango <- log(range(datos[datos > umbral], na.rm=T), base)
  escala <- seq(from=rango[1], to=rango[2], length.out=nIntervalos-1)
  escala <- unique(c(min(datos, na.rm=T), round(base^escala, nDigitos)))
  
  return (crearEscala(escala=escala, colores=colores, brewerPal = brewerPal, invertirPaleta=invertirPaleta, continuo=continuo))
}

crearEscalaEnQuantiles <- function(datos, nDigitos=1, nIntervalos=9, probs=seq(from = 0, to=1, length.out = nIntervalos+1), 
                                   colores=NULL, brewerPal='Spectral', invertirPaleta=brewerPal %in% paletasInvertidas, continuo=F) {
  nIntervalos <- length(probs)-1
  escalaOrig <- quantile(datos, probs, na.rm = T)
  
  escala <- escalaOrig
  escala[1] <- achicarToNDigitos(escala[1], nDigitos = nDigitos)
  escala[length(escala)] <- agrandarToNDigitos(escala[length(escala)], nDigitos = nDigitos)
  escala <- round(escala, nDigitos)
  
  nDigitosOrig <- nDigitos
  maxNDigitos <- nDigitosOrig * 3
  dups <- duplicated(escala)
  while (any(dups) && nDigitos <= maxNDigitos) {
    escala <- escalaOrig
    nDigitos <- nDigitos + 1
    escala[1] <- achicarToNDigitos(escala[1], nDigitos = nDigitos)
    escala[length(escala)] <- agrandarToNDigitos(escala[length(escala)], nDigitos = nDigitos)
    escala <- round(escala, nDigitos)
    dups <- duplicated(escala)
  }
  
  if (any(dups)) escala <- escala[!dups]
  escala <- as.numeric(escala)
  
  #i <- 2
  #iUltimoNoDuplicado <- 1
  #while (any(dups) && i < length(dups)) {
    # busco el último duplicado desde donde estoy
  #  while (i < length(dups) & dups[i])
  #    i <- i + 1
    
  #  probs <- c(probs[1:iUltimoNoDuplicado], seq(from=probs[i], to=1, length.out = nIntervalos - iUltimoNoDuplicado))
  #  escala <- quantile(datos, probs, na.rm = T)
  #  escala[1] <- achicarToNDigitos(escala[1], nDigitos = nDigitos)
  #  escala[length(escala)] <- agrandarToNDigitos(escala[length(escala)], nDigitos = nDigitos)
  #  escala <- round(escala, nDigitos)
  #  dups <- duplicated(escala)
    
  #  while (i < length(dups) & !dups[i])
  #    i <- i + 1
    
  #  iUltimoNoDuplicado <- i - 1
  #}
  #escala <- unique(quantile(datos, probs, na.rm = T))
  #escala[1] <- achicarToNDigitos(escala[1], nDigitos = nDigitos)
  #escala[length(escala)] <- agrandarToNDigitos(escala[length(escala)], nDigitos = nDigitos)
  #escala <- round(escala, nDigitos)
  
  return(crearEscala(escala=escala, colores=colores, brewerPal = brewerPal, invertirPaleta=invertirPaleta, continuo=continuo))
}

crearEscalaEquiespaciadaEnQuantiles <- function(datos, nDigitos=1, nIntervalos=9, colores=NULL, brewerPal='Spectral', invertirPaleta=brewerPal %in% paletasInvertidas, continuo=F) {
  return(crearEscalaEnQuantiles(datos = datos, nDigitos = nDigitos, nIntervalos = nIntervalos, colores = colores, brewerPal = brewerPal, invertirPaleta=invertirPaleta, continuo = continuo))
}

crearEscalaEquiespaciadaEnQuintilesDistinguirMenoresOIgualesAUmbral <- function(datos, umbral, nDigitos=1, nIntervalos=9, colores=NULL, brewerPal='Spectral', invertirPaleta=brewerPal %in% paletasInvertidas, continuo=F) {
  # Crea una escala que tiene en el primer intervalo todos los valores menores o iguales a
  # umbral, y en el resto de los intervalos es equiespaciada desde el mínimo valor en datos mayor a umbral 
  # hasta el máximo valor
  datosNoNA <- datos[!is.na(datos)]
  datosNoNAMayoresAUmbral <- datosNoNA[datosNoNA > umbral]
  
  probs <- seq(from = 0, to=1, length.out = nIntervalos - 1)
  escala <- unique(round(c(min(datosNoNA), quantile(datosNoNAMayoresAUmbral, probs)), nDigitos))  
  
  return (crearEscala(escala=escala, colores=colores, brewerPal = brewerPal, invertirPaleta=invertirPaleta, continuo=continuo))
}

crearEscalaEnBaseASubconjuntoDeLosDatos <- function(datos, iSubConjunto, nDigitos=1, nIntervalos=9, colores=NULL, brewerPal='Spectral', invertirPaleta=brewerPal %in% paletasInvertidas, continuo=F) {
  rangoSubConjunto <- range(datos[iSubConjunto], na.rm=T)
  rango <- range(datos, na.rm=T)
  
  escala <- unique(round(c(achicarToNDigitos(rango[1], nDigitos = nDigitos),
                   seq(from=achicarToNDigitos(rangoSubConjunto[1]), to=agrandarToNDigitos(rangoSubConjunto[2]), length.out=nIntervalos-1),
                   agrandarToNDigitos(rango[2], nDigitos = nDigitos)
                   ), nDigitos))
  return (crearEscala(escala=escala, colores=colores, brewerPal = brewerPal, invertirPaleta=invertirPaleta, continuo=continuo))
}

ajustarExtremosEscala <- function(escala, datos, nDigitos=1, redondear=TRUE) {
  if (length(escala$escala) > 1) {
    rango <- rangoExtendidoANDigitos(x = datos, nDigitos = nDigitos, na.rm = T)
    nEscala <- length(escala$escala)
    if (escala$escala[1] > rango[1]) escala$escala[1] <- rango[1]
    if ((escala$colores[length(escala$colores)] == '' | escala$continuo) && escala$escala[nEscala] < rango[2]) escala$escala[nEscala] <- rango[2]
    if (redondear) escala$escala <- c(achicarToNDigitos(escala$escala[1], nDigitos = nDigitos), 
                                      round(escala$escala[2:(nEscala-1)], digits = nDigitos), 
                                      agrandarToNDigitos(escala$escala[nEscala], nDigitos = nDigitos))
    
    iUnicos <- !duplicated(escala$escala)
    escala$escala <- escala$escala[iUnicos]
    escala$colores <- escala$colores[iUnicos]
  }
  return(escala)
}

getXYLims <- function(spObjs, resXImagenes=640, resYImagenes=NULL, ejesXYLatLong=TRUE, factorMargen = 0.025) {
  # obtengo el rectángulo que encierra a la grilla u observaciones y dejo un borde de 2.5% para cada lado
  p4str <- proj4string(spObjs[[1]])
  xLim <- spObjs[[1]]@bbox[1,]
  yLim <- spObjs[[1]]@bbox[2,]

  if (length(spObjs) > 1) {
    for (i in 2:length(spObjs)) {
      if (!is.null(spObjs[[i]])) {
        xLimI <- spObjs[[i]]@bbox[1,]
        yLimI <- spObjs[[i]]@bbox[2,]
        
        if (xLim[1] > xLimI[1]) xLim[1] <- xLimI[1]
        if (xLim[2] < xLimI[2]) xLim[2] <- xLimI[2]
        if (yLim[1] > yLimI[1]) yLim[1] <- yLimI[1]
        if (yLim[2] < yLimI[2]) yLim[2] <- yLimI[2]
      }
    }
  }
  
  if (factorMargen > 1E-3) {
    deltaX <- (xLim[2] - xLim[1]) * factorMargen
    xLim[1] <- trunc(xLim[1] - deltaX)
    xLim[2] <- ceiling(xLim[2] + deltaX)
    deltaY <- (yLim[2] - yLim[1]) * factorMargen
    yLim[1] <- trunc(yLim[1] - deltaY)
    yLim[2] <- ceiling(yLim[2] + deltaY)
  }
  
  if (is.null(resYImagenes) || resYImagenes <= 0) { resYImagenes <- round(resXImagenes * (yLim[2] - yLim[1]) / (xLim[2] - xLim[1]))
  } else { if (is.null(resXImagenes) || resXImagenes <= 0) { resXImagenes <- round(resYImagenes * (xLim[2] - xLim[1]) / (yLim[2] - yLim[1])) } }  

  if (ejesXYLatLong) {
    # Hago el "cuadrado" en lat/long
    p4strLatLong <- '+proj=longlat +datum=WGS84'
    coordsBB <- matrix(data=c(xLim[1], yLim[1], #abajo, izquierda
                              xLim[1], yLim[2], #arriba, izquierda
                              xLim[2], yLim[2], #arriba, derecha
                              xLim[2], yLim[1], #abajo, derecha
                              xLim[1], yLim[1]), #abajo, izquierda
                       ncol=2, byrow=T)
    p <- Polygon(coordsBB)
    ps <- Polygons(list(p), 1)
    sps <- SpatialPolygons(list(ps), proj4string=CRS(p4str))
    spsLatLong <- spTransform(sps, CRS(p4strLatLong))
    
    # Redondeo para sacar los decimales
    lonMin <- floor(spsLatLong@bbox[1, 1])
    lonMax <- ceiling(spsLatLong@bbox[1, 2])
    latMin <- floor(spsLatLong@bbox[2, 1])
    latMax <- ceiling(spsLatLong@bbox[2, 2])
    
    # Creo las rectas en Long/Lat verticales y horizontales en cada grado entero
    xLineas <- seq(from=lonMin, to=lonMax, by=1)
    yLineas <- seq(from=latMin, to=latMax, by=0.2)
    
    lineasVerticales <- list()
    length(lineasVerticales) <- length(xLineas)
    i <- 1
    for (lon in xLineas) {
      lineaVertical <- Line(expand.grid(x=lon, y=yLineas))
      lineaVertical <- Lines(list(lineaVertical), ID=paste('x',as.character(lon),sep=''))
      lineasVerticales[[i]] <- lineaVertical
      i <- i + 1
    }
    
    xLineas <- seq(from=lonMin, to=lonMax, by=0.2)
    yLineas <- seq(from=latMin, to=latMax, by=1)
    
    lineasHorizontales <- list()
    length(lineasHorizontales) <- length(yLineas)
    i <- 1
    for (lat in yLineas) {
      lineaHorizontal <- Line(expand.grid(x=xLineas, y=lat))
      lineaHorizontal <- Lines(list(lineaHorizontal), ID=paste('y',as.character(lat),sep=''))
      lineasHorizontales[[i]] <- lineaHorizontal
      i <- i + 1
    }
    lineasLatLong <- SpatialLines(c(lineasVerticales, lineasHorizontales), proj4string=CRS(p4strLatLong))
    lineasLatLong <- spTransform(lineasLatLong, CRS(p4str))
    # Intersecto las lineas horizontales y verticales con el cuadrado lat/long proyectado
    # para quedarme solo con la parte dentro del área del gráfico
    lineasLatLong <- gIntersection(lineasLatLong, sps, byid=T)
  
    # Obtengo la intersección de cada recta con el marco del área de ploteo del mapa
    # y el texto a mostrar en cada intersección, W o E para Lon y S o N para Lat y el grado entero
    # Borde inferior del marco de ploteo
    linea <- SpatialLines(list(Lines(list(Line(expand.grid(x=xLim, y=yLim[1]))), ID='aux')), proj4string=CRS(p4str))
    cortesEjeX <- gIntersection(lineasLatLong, linea)
    breaksEjeX <- cortesEjeX@coords[,1]
    row.names(cortesEjeX) <- 1:length(cortesEjeX)
    cortesEjeX <- spTransform(cortesEjeX, CRS(p4strLatLong))
    
    textoEjeX <- character(0)
    lons <- cortesEjeX@coords[,1]
    i<-1
    for (i in 1:length(lons)) {
      if (lons[i] < 0) { textoEjeX <- c(textoEjeX, sprintf('%dW', -as.integer(lons[i])))
      } else { textoEjeX <- c(textoEjeX, sprintf('%dE', as.integer(lons[i]))) }
    }    
    
    # Borde izquierdo del marco de ploteo
    linea <- SpatialLines(list(Lines(list(Line(expand.grid(x=xLim[1], y=yLim))), ID='aux')), proj4string=CRS(p4str))
    cortesEjeY <- gIntersection(lineasLatLong, linea)
    breaksEjeY <- cortesEjeY@coords[,2]
    row.names(cortesEjeY) <- 1:length(cortesEjeY)
    cortesEjeY <- spTransform(cortesEjeY, CRS(p4strLatLong))
    textoEjeY <- character(0)
    lats <- cortesEjeY@coords[,2]
    for (i in 1:length(lats)) {
      if (lats[i] < 0) { textoEjeY <- c(textoEjeY, sprintf('%dS', -as.integer(lats[i])))
      } else { textoEjeY <- c(textoEjeY, sprintf('%dN', as.integer(lats[i]))) }
    }    

    lineasLatLong <- SpatialLinesDataFrame(sl=lineasLatLong, data=data.frame(dummy=rep(1, length(lineasLatLong))), match.ID=F)
    ejesLatLong <- list(lineasLatLong=lineasLatLong, breaksEjeX=breaksEjeX, textoEjeX=textoEjeX, 
                        breaksEjeY=breaksEjeY, textoEjeY=textoEjeY)
  } else {
    ejesLatLong <- NULL
  }
  return(list(xLim=xLim, yLim=yLim, resXImagenes=resXImagenes, resYImagenes=resYImagenes, ejesLatLong=ejesLatLong))
}

getLayoutSHP <- function(shp) {
  return(list("sp.polygons", shp, zcol=1, fill=NA, col=rgb(96, 96, 96, maxColorValue=255), first=F))
}

aplicarOpcionesAMapa <- function(p, xyLims, shpBase, dibujarEscala=T, dibujarEjes=T, tamanioFuenteEjes=15, 
                                 tamanioFuenteTitulo=14, tamanioFuenteSubtitulo=12,
                                 titulo='', subtitulo='', colorFillSHPBase=NA, escalaGraficos=1, 
                                 puntosAResaltar=NULL, colorResalto=rgb(255, 0, 0, maxColorValue=255),
                                 tamanioResalto=0.8, widthPx) {
  if (!is.null(shpBase)) {
    if (!rgeos::gIsValid(shpBase)) shpBase <- gBuffer(shpBase, byid=TRUE, width=0)
    if ('SpatialPolygonsDataFrame' %in% class(shpBase)) {
      shpBase@data$id <- rownames(shpBase@data)  
    } else {
      shpBase <- SpatialPolygonsDataFrame(Sr = shpBase, data = data.frame(id=1:length(shpBase)))
    }
    
    shpF <- fortify(shpBase, region="id")
    if (is.na(colorFillSHPBase)) {
      p <- p + geom_path(data=shpF, mapping=aes(x=long, y=lat, group=group, z=NULL), color=rgb(96, 96, 96, maxColorValue=255))
    } else {
      p <- p + geom_polygon(data=shpF, mapping=aes(x=long, y=lat, group=group, z=NULL), color=rgb(96, 96, 96, maxColorValue=255), fill=colorFillSHPBase,
                            show.legend = FALSE)
    }
  }
  
  if (!is.null(puntosAResaltar)) {
    minLargo <- min(diff(xyLims$xLim), diff(xyLims$yLim))
    poligonosResalto <- gBuffer(spgeom = puntosAResaltar, byid = T, width = minLargo * 0.05)
    # poligonosResalto <- gBuffer(spgeom = puntosAResaltar, byid = T, width = maxDist, quadsegs = 16)
    if ('SpatialPolygonsDataFrame' %in% class(poligonosResalto)) { poligonosResalto@data$id <- rownames(poligonosResalto@data)
    } else { poligonosResalto <- SpatialPolygonsDataFrame(Sr = poligonosResalto, data = data.frame(id=1:length(poligonosResalto))) }
    shpF <- fortify(poligonosResalto, region="id")
    p <- p + geom_polygon(data=shpF, mapping=aes(x=long, y=lat, group=group, z=NULL), color=colorResalto, fill=NA, show.legend = FALSE, size=tamanioResalto)
  }
  
  if (!is.null(xyLims$ejesLatLong)) {
    xyLims$ejesLatLong$lineasLatLong@data$id <- rownames(xyLims$ejesLatLong$lineasLatLong@data)
    lineasLatLongF <- fortify(xyLims$ejesLatLong$lineasLatLong, region="id")
    
    p <- p + 
      geom_line(data=lineasLatLongF, mapping=aes(x=long, y=lat, group=group, z=NULL), color=gray(level=0.2, alpha=0.4)) +
      scale_x_continuous(breaks=xyLims$ejesLatLong$breaksEjeX, labels=xyLims$ejesLatLong$textoEjeX) +
      scale_y_continuous(breaks=xyLims$ejesLatLong$breaksEjeY, labels=xyLims$ejesLatLong$textoEjeY) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  }
  
  if (!dibujarEjes) {
    p <- p + theme(axis.line=element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank(),
                   axis.ticks=element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank())  
  }
  
  p <- p + theme(text = element_text(size=tamanioFuenteEjes*escalaGraficos))
  
  # dibujarEscala <- T
  if (!dibujarEscala) p <- p + theme(legend.position="none")
  
  wrap_strings <- function(vector_of_strings, width) { 
    vector_of_strings <- strsplit(vector_of_strings, split = '\n')
    vector_of_strings <- sapply(vector_of_strings, FUN=function(x) {paste(strwrap(x, width=width), collapse="\n")}, USE.NAMES = FALSE)    
    return(paste(vector_of_strings, collapse = '\n'))
  } 
  
  maxWidthChars <- trunc(65 * widthPx / 630)
  if (titulo != '') {
    p <- p + labs(title = wrap_strings(titulo, maxWidthChars)) + 
          theme(plot.title = element_text(size = tamanioFuenteTitulo * escalaGraficos, face = "bold", colour = "black", hjust=0.5))
  }
  
  if (subtitulo != '') {  
    p <- p + labs(subtitle = wrap_strings(subtitulo, maxWidthChars)) + 
          theme(plot.subtitle = element_text(size = tamanioFuenteSubtitulo * escalaGraficos, face="italic", colour = "black", hjust=0.5))
  }
  
  return(p)
}

mapearPuntos <- function(observaciones, layoutSHP=NULL, escala, nomArchResultados, xyLims=NULL,
                         dibujarEscala=T, dibujarEjes=T, zcol=1, titulo='') {
  # TO-DO: agregar título al gráfico
  oldSciPen <- getOption("scipen")
  options(scipen=15)
  if (is.null(xyLims)) xyLims <- getXYLims(c(observaciones))
  if (dibujarEjes) { colEjes <- 'black' } else { colEjes <- 'transparent' }
  parSettings <- list(axis.line=list(col=colEjes))
  
  # bubble(obj=observaciones[!is.na(observaciones$value),], zcol='value', scales=list(draw=dibujarEjes), sp.layout=layoutSHP, xlim=xLim, ylim=yLim, do.sqrt=T, col=escalaObservaciones$colores[1:nrow(escalaObservaciones)-1])
  if (!dibujarEscala) {
    plot <- spplot(observaciones, zcol=zcol, col.regions=escala$colores[1:(length(escala$escala)-1)], 
                   cuts=escala$escala, scales=list(draw=dibujarEjes), sp.layout=layoutSHP, 
                   xlim=xyLims$xLim, ylim=xyLims$yLim, colorkey=FALSE, key=list(lines=TRUE, col="transparent"), par.settings=parSettings)
  } else {
    plot <- spplot(observaciones, 'value', col.regions=escala$colores[1:(length(escala$escala)-1)],
                   cuts=escala$escala, scales=list(draw=dibujarEjes), sp.layout=layoutSHP, 
                   xlim=xyLims$xLim, ylim=xyLims$yLim, colorkey=TRUE, par.settings=parSettings)
  }
  png(nomArchResultados, height=xyLims$resYImagenes, width=xyLims$resXImagenes, type='cairo')
  tryCatch(expr = print(plot), finally = dev.off())
  if (interactive()) { print(plot) }

  options(scipen=oldSciPen)
}

mapearGrilla <- function(grilla, layoutSHP=NULL, escala, nomArchResultados, xyLims=NULL, 
                         dibujarEscala=T, dibujarEjes=T, zcol=1, contour=F, labels=F) {
  oldSciPen <- getOption("scipen")
  options(scipen=15)
  if (is.null(xyLims)) xyLims <- getXYLims(c(grilla))
  if (dibujarEjes) { colEjes <- 'black' } else { colEjes <- 'transparent' }
  parSettings <- list(axis.line=list(col=colEjes))
  
  # contour=T, labels=T, para dibujar las líneas de contorno y etiquetarlas con su valor  
  plotPred <- spplot(obj=grilla, zcol=zcol, col.regions=escala$colores[1:(length(escala$escala)-1)], 
                     at=escala$escala, pretty=F, scales=list(draw=dibujarEjes), sp.layout=layoutSHP, 
                     xlim=xyLims$xLim, ylim=xyLims$yLim, colorkey=dibujarEscala, par.settings=parSettings,
                     contour=contour, labels=labels)
  png(nomArchResultados, height=xyLims$resYImagenes, width=xyLims$resXImagenes, type='cairo')
  tryCatch(expr = print(plotPred), finally = dev.off())
  if (interactive()) { print(plotPred) }
  
  options(scipen=oldSciPen)
}

mapearPuntosConEtiquetasGGPlot <- function(puntos, shpBase=NULL, nomArchResultados=NULL, xyLims=NULL, 
                                           dibujarEjes=T, zcol=1, DPI=90, widthPx=630, heightPx=630,
                                           tamanioFuentePuntos=3, nDigitos=2, dibujarPuntos=T, 
                                           coloresPuntos="black", coloresTexto="black", tamaniosPuntos=3, 
                                           escala=NULL, dibujar=interactive(), titulo='', subtitulo='') {
  if (exists(".Random.seed", .GlobalEnv)) { set.seed(31) }
  oldseed <- .GlobalEnv$.Random.seed
  set.seed(31)
  
  oldSciPen <- getOption("scipen")
  options(scipen=15)
  if (is.null(xyLims)) xyLims <- getXYLims(spObjs = c(puntos, shpBase))
  if (dibujarPuntos) { justV <- -1
  } else { justV <- 0 }
  
  iNoNA <- !is.na(puntos@data[,zcol])
  todosNA <- !any(iNoNA)
  v <- puntos@data[iNoNA, zcol]
  
  if (!is.null(escala) && is.numeric(v)) {
    nEscala <- length(escala$escala)
    if (nEscala == 1) {
      v <- factor(v)
      valMax <- escala$escala[nEscala]
      ultimoI <- 1
      continuo <- FALSE
    } else {
      # TO-DO: escalas de colores continuas todavía no están implenentadas
      continuo <- F
      if (continuo || escala$colores[nEscala] == "") { 
        valMax <- escala$escala[nEscala]
        ultimoI <- nEscala - 1
      } else { 
        valMax <- Inf 
        ultimoI <- nEscala
      }
      
      v <- cut(v, breaks=c(escala$escala[1:ultimoI], valMax), right=FALSE, include.lowest=!is.infinite(valMax))
    }
    #coloresTexto <- escala$colores[factor(v)]
    coloresPuntos <- escala$colores[factor(v)]
  }
  
  # TO-DO: Revisar esto. Puede ser mejor usar un color distinto para los puntos sin datos
  # puntos <- puntos[!is.na(puntos@data[,zcol]), ]

  coords <- coordinates(puntos)[iNoNA, , drop=F]
  value <- puntos@data[iNoNA, zcol]
  if (is.double(value) && !is.POSIXt(value)) { value <- round(x=value, nDigitos) }
  
  df <- data.frame(x=coords[,1], y=coords[,2], value=value)
  
  p <- ggplot(aes(x=x, y=y, z=value), data=df) + 
       coord_equal(ratio=1, xlim=xyLims$xLim, ylim=xyLims$yLim) +
       #labs(x = "Este[m]", y = "Norte[m]", fill = "") + 
       labs(x = "", y = "", fill = "") + 
       theme(panel.background=element_blank())

  # Dibujar escala no se usa, se pasa con true para que no haga nada
  p <- aplicarOpcionesAMapa(p=p, xyLims=xyLims, shpBase=shpBase, dibujarEscala=T, dibujarEjes=dibujarEjes, 
                            titulo=titulo, subtitulo=subtitulo, widthPx = widthPx)  
  
  if (dibujarPuntos) p <- p + geom_point(data=df, aes(x=x, y=y), colour=coloresPuntos, size=tamaniosPuntos)
  #p <- p + geom_text(data=df, aes(label=value), vjust=justV, size=tamanioFuentePuntos, colour=coloresTexto)
  p <- p + geom_text_repel(data=df, aes(label=value), size=tamanioFuentePuntos, colour=coloresTexto)
  
  #if (!is.null(escala) && is.numeric(value)) {
  #  breaks <- levels(v)
  #  p <- p + scale_fill_manual(breaks=rev(breaks), drop=F, 
  #                             labels=rev(breaks),
  #                             values=c(escala$colores[1:ultimoI]))
  #}
  
  if (dibujar) print(p)
  if (!is.null(nomArchResultados)) {
    path <- dirname(nomArchResultados)
    if (!file.exists(path)) dir.create(path, showWarnings=F, recursive=T)
    ggsave(p, file=nomArchResultados, dpi=DPI, width = widthPx / DPI, height = heightPx / DPI, units = 'in', type='cairo')
  }
  options(scipen=oldSciPen)
  
  
  .GlobalEnv$.Random.seed <- oldseed 
  
  return(p)
}

mapearPuntosGGPlot <- function(puntos, shpBase=NULL, nomArchResultados=NULL, xyLims=NULL, dibujarEjes=T, zcol=1, DPI=90, widthPx=630, heightPx=630, 
                               tamaniosPuntos=5, dibujarTexto=F, tamanioFuentePuntos=3, tamanioFuenteEjes=15, nDigitos=1, escala=NULL, 
                               dibujar=interactive(), titulo='', subtitulo='', colorFillSHPBase=NA, contornearPuntos=FALSE, continuo, 
                               alturaEscalaContinua=unit(0.1, 'npc'), escalaGraficos = 1, puntosAResaltar=NULL) {
  oldSciPen <- getOption("scipen")
  options(scipen=15)
  if (is.null(xyLims)) {
    if (is.null(shpBase)) { xyLims <- getXYLims(spObjs = list(puntos), ejesXYLatLong = F)
    } else { xyLims <- getXYLims(spObjs = list(puntos, shpBase), ejesXYLatLong = F) }
  }
  iNoNA <- !is.na(puntos@data[,zcol])
  todosNA <- !any(iNoNA)
  v <- puntos@data[iNoNA,zcol]
  
  # continuo <- T
  if (missing(continuo)) continuo <- !is.null(escala) && !is.null(escala$continuo) && escala$continuo
  
  if (is.null(escala)) { 
    if (!todosNA) { escala <- crearEscalaEquiespaciada(datos=v, continuo = continuo)
    } else { escala <- crearEscalaEquiespaciada(datos=c(0, 1), continuo = continuo) }
  }
  
  nEscala <- length(escala$escala)
  if (nEscala == 1) {
    v <- factor(v)
    valMax <- escala$escala[nEscala]
    ultimoI <- 1
    continuo <- FALSE
  } else {
    if (continuo || escala$colores[nEscala] == "") { 
      valMax <- escala$escala[nEscala]
      ultimoI <- nEscala - 1
    } else { 
      valMax <- Inf 
      ultimoI <- nEscala
    }
    
    if (!continuo) { v <- cut(v, breaks=c(escala$escala[1:ultimoI], valMax), right=FALSE, include.lowest=!is.infinite(valMax), dig.lab=15, ordered_result=T) }
  }
    
  coords <- coordinates(puntos)[iNoNA,,drop=F]
  if (length(tamaniosPuntos) > 1) tamaniosPuntos <- tamaniosPuntos[iNoNA]
  df <- data.frame(x=as.numeric(coords[,1]), y=as.numeric(coords[,2]), value = v, size = tamaniosPuntos * escalaGraficos)
  # df$size <- tamaniosPuntos * escalaGraficos

  p <- ggplot(aes(x=x, y=y), data=df) + 
              labs(x = '', y = '', colour = '', size='') + 
              coord_equal(ratio=1, xlim=xyLims$xLim, ylim=xyLims$yLim) +
              #labs(x = "Este[m]", y = "Norte[m]", fill = "") + 
              theme(panel.background=element_blank())
  
  if (continuo && nEscala > 1) {
    vals <- ggplot2:::rescale01(escala$escala)
    alturaEscalaContinuaEscalada <- unit(x=as.numeric(alturaEscalaContinua) * escalaGraficos, units = attr(alturaEscalaContinua, 'unit'))
    p <- p + scale_color_gradientn(colours=escala$colores, values=vals, limits=c(escala$escala[1], valMax),
                                   breaks=escala$escala, na.value="gray95") + theme(legend.key.height=alturaEscalaContinuaEscalada)
  } else {
    breaks <- rev(levels(v))
    if (length(breaks) > 11) { labels <- rev(escala$escala[1:ultimoI])
    } else { labels <- breaks }
    labels <- gsub(pattern = ',', replacement = ', ', x = labels, fixed = T)
    #breaks <- paste('[', escala$escala[1:ultimoI], ',', c(escala$escala[2:ultimoI], as.character(valMax)), ')', sep='')
    p <- p + scale_colour_manual(breaks=breaks, drop=F, labels=labels, values=escala$colores[1:ultimoI], na.value="gray95")
  }

  # Dibujar escala no se usa, se pasa con true para que no haga nada
  #dibujarEjes <- T
  #subtitulo <- ''
  p <- aplicarOpcionesAMapa(p=p, xyLims=xyLims, shpBase=shpBase, dibujarEscala=T, dibujarEjes=dibujarEjes, 
                            tamanioFuenteEjes=tamanioFuenteEjes * escalaGraficos, titulo=titulo, subtitulo=subtitulo,
                            colorFillSHPBase=colorFillSHPBase, puntosAResaltar = puntosAResaltar, widthPx = widthPx)
  if (!todosNA) {
    if (contornearPuntos)
      p <- p + geom_point(colour='black', size=tamaniosPuntos * escalaGraficos + 1, show.legend = F)
    
    if (length(tamaniosPuntos) == 1) {
      p <- p + geom_point(aes(colour=value), size=tamaniosPuntos * escalaGraficos)
    } else {
      p <- p + 
        geom_point(aes(colour=value, size=factor(size))) +
        scale_size_manual(values=sort(unique(tamaniosPuntos)) * escalaGraficos, guide='none')
    }
    
    if (dibujarTexto) {
      labels <- round(puntos@data[iNoNA, zcol], nDigitos)
      p <- p + geom_text_repel(data=df, aes(label=labels), size=tamanioFuentePuntos * escalaGraficos)
      #p <- p + geom_text(data=df, aes(label=labels), vjust=-1, size=tamanioFuentePuntos * escalaGraficos)
    }
  }

  if (dibujar) print(p)
  if (!is.null(nomArchResultados)) {
    path <- dirname(nomArchResultados)
    if (!file.exists(path)) dir.create(path, showWarnings=F, recursive=T)
    ggsave(p, file=nomArchResultados, dpi=DPI * escalaGraficos, width = (widthPx / DPI) * escalaGraficos, 
           height = (heightPx / DPI) * escalaGraficos, units = 'in', type='cairo')
  }
  options(scipen=oldSciPen)
  
  return(p)
}

map_aspect = function(x, y) {
  x.center <- sum(range(x)) / 2
  y.center <- sum(range(y)) / 2
  x.dist <- ggplot2:::dist_central_angle(x.center + c(-0.5, 0.5), rep(y.center, 2))
  y.dist <- ggplot2:::dist_central_angle(rep(x.center, 2), y.center + c(-0.5, 0.5))
  y.dist / x.dist
}

mapearGrillaGGPlot <- function(grilla, shpBase=NULL, escala=NULL, nomArchResultados=NULL, xyLims=NULL,
                               dibujarEscala=TRUE, dibujarEjes=TRUE, zcol=1, isolineas=FALSE, DPI=90,
                               widthPx=630, heightPx=630, dibujar=interactive(), titulo='',
                               subtitulo = '', continuo, alturaEscalaContinua=unit(0.1, 'npc'),
                               dibujarPuntosObservaciones=FALSE, coordsObservaciones=NULL, 
                               tamaniosPuntos = 0.8, tamanioFuentePuntos = 3, puntosAResaltar=NULL) {
  #grilla <- coarsenGrid(grilla, coarse = 6)
  if (is.null(xyLims)) {
    if (is.null(shpBase)) { xyLims <- getXYLims(spObjs = list(grilla), ejesXYLatLong = F)
    } else { xyLims <- getXYLims(spObjs = list(grilla, shpBase), ejesXYLatLong = F, factorMargen = 0) }
  }

  # continuo <- T
  if (missing(continuo)) continuo <- !is.null(escala) && !is.null(escala$continuo) && escala$continuo
    
  iNoNa <- !is.na(grilla@data[, zcol])
  todosNA <- !any(iNoNa)
  v <- grilla@data[iNoNa, zcol]
  if (is.null(escala)) { 
    if (!todosNA) { escala <- crearEscalaEquiespaciada(datos=v, continuo = continuo)
    } else { escala <- crearEscalaEquiespaciada(datos=c(0, 1), continuo = continuo) }
  }

  nEscala <- length(escala$escala)
  if (nEscala == 1) {
    valMax <- escala$escala[1] + 1
    v <- factor(v)
    ultimoI <- 1
    continuo <- FALSE
  } else {
    if (continuo || escala$colores[nEscala] == "") { 
      valMax <- escala$escala[nEscala]
      ultimoI <- nEscala - 1
    } else { 
      valMax <- Inf 
      ultimoI <- nEscala
    }
  
    if (!continuo) { v <- cut(v, breaks=c(escala$escala[1:ultimoI], valMax), right=FALSE, include.lowest=!is.infinite(valMax), dig.lab=15, ordered_result=T) }
  }
  
  coords <- coordinates(grilla)[iNoNa, , drop=F]
  df <- data.frame(x=coords[, 1], y=coords[, 2], value=v)
  if (is.projected(grilla)) { ratio <- 1 
  } else { ratio <- map_aspect(df$x, df$y) }
  
  p <- ggplot(aes(x=x, y=y, z=value), data=df) +
       coord_equal(ratio=ratio, xlim=xyLims$xLim, ylim=xyLims$yLim) +
       # labs(x = "Este[m]", y = "Norte[m]", fill = "") +
       labs(x = "", y = "", fill = "") +
       theme(panel.background=element_blank())

  if (gridded(grilla)) {
    # p <- p + geom_tile(data=df, aes(fill=value), show.legend = TRUE)
    p <- p + geom_raster(data=df, aes(fill=value), show.legend = TRUE)
  } else { 
    #p <- p + geom_tile(data=df, aes(fill=value), show.legend = TRUE)
    p <- p + geom_point(data=df, aes(color=value), show.legend = TRUE)
    # stop('mapearGrillaGGPlot no implementado para gridded(grilla) == F') 
  }
  
  if (!todosNA){
    if (continuo && isolineas && nrow(df) > 0 && var(grilla@data[iNoNa,zcol]) > 1E-6) {
      if (!is.null(escala$iniciosIntervalosIsoLineas)) { cortesIsolineas <- escala$iniciosIntervalosIsoLineas
      } else  { 
        nPuntos <- 7
        cortesIsolineas <- round(quantile(df$value, seq(from = 1/nPuntos, to = 1 -1/nPuntos, length.out = nPuntos - 2)), 1)
      }
        
      p <- p + geom_contour(aes(colour = ..level..), breaks=cortesIsolineas, color='gray30', na.rm=T, show.legend=T)# + scale_colour_gradient(guide = 'none')
      instant_pkgs("directlabels")
      #p <- direct.label(p, list("far.from.others.borders", "calc.boxes", "enlarge.box", 
      #                  hjust = 1, vjust = 1, box.color = NA, fill = "transparent", "draw.rects"))
      p <- direct.label(p, list("bottom.pieces", colour='gray30', size=3))
    }
  }
  
  if (continuo) {
    vals <- ggplot2:::rescale01(escala$escala)
    p <- p + scale_fill_gradientn(colours=escala$colores, values=vals, limits=c(escala$escala[1], valMax),
                                  breaks=escala$escala, na.value="gray95") + theme(legend.key.height=alturaEscalaContinua)
  } else {
    breaks <- rev(levels(v))
    #if (length(breaks) > 18) { labels <- rev(escala$escala[1:ultimoI])
    #} else { 
      labels <- breaks 
    #}
    labels <- gsub(pattern = ',', replacement = ', ', x = labels, fixed = T)
    #breaks <- paste('[', escala$escala[1:ultimoI], ',', c(escala$escala[2:ultimoI], as.character(valMax)), ')', sep='')
    p <- p + scale_fill_manual(breaks=breaks, drop=F, labels=labels, values=escala$colores[1:ultimoI], na.value="gray95")
  }    
  
  p <- aplicarOpcionesAMapa(p=p, xyLims=xyLims, shpBase=shpBase, dibujarEscala=dibujarEscala, dibujarEjes=dibujarEjes, 
                            titulo=titulo, subtitulo=subtitulo, puntosAResaltar = puntosAResaltar, widthPx = widthPx)
  
  if (dibujarPuntosObservaciones & !is.null(coordsObservaciones)) {
    zColObs <- max(which(colnames(coordsObservaciones@data) == 'value'), 1)
    if (!is.null(shpBase)) { coordsObsSobreShpBase <- coordsObservaciones[!is.na(over(coordsObservaciones, geometry(shpBase))),]
    } else { coordsObsSobreShpBase <- coordsObservaciones }
    coordsObsSobreShpBase <- coordsObsSobreShpBase[!is.na(coordsObsSobreShpBase@data[, zColObs]), ]
    
    if (length(coordsObsSobreShpBase) > 60) {
      iMuestras <- muestrearEnCuadrantesYECDF(sp = coordsObsSobreShpBase, size = 60, nCuadrantesX = 2, nCuadrantesZ = 10, zcol = zcol)
      coordsObsSobreShpBase <- coordsObsSobreShpBase[iMuestras,]
    }
    coordsPuntos <- coordinates(coordsObsSobreShpBase)
    dfPuntos <- data.frame(x=coordsPuntos[, 1], y=coordsPuntos[, 2], value=round(coordsObsSobreShpBase@data[, zColObs], 1))
    tamaniosPuntos = 1
    p <- p + geom_point(data=dfPuntos, aes(x=x, y=y), colour='black', size=tamaniosPuntos) + 
         geom_text_repel(data=dfPuntos, aes(label=value), size=tamanioFuentePuntos, colour='black')
  }
  
  oldSciPen <- getOption("scipen")
  options(scipen=15)
  if (dibujar) print(p)
  if (!is.null(nomArchResultados)) {
    path <- dirname(nomArchResultados)
    if (!file.exists(path)) dir.create(path, showWarnings=F, recursive=T)
    
    ggsave(p, file=nomArchResultados, dpi=DPI, width = widthPx / DPI, height = heightPx / DPI, units = 'in', type='cairo')
  }
  options(scipen=oldSciPen)
  
  return(p)
}

mapearEx <- function(observaciones, grilla, shpMask=NULL, escalaObservaciones, escalaSD,
                     nomArchResultados, resXImagenes=640, resYImagenes=NULL, 
                     graficarObservaciones=T, graficarPrediccion=T, graficarVarianza=F, graficarSD=F, 
                     dibujarEscala=T, dibujarEjes=T) {
  spObjs <- list()
  if (!is.null(shpMask)) spObjs[[length(spObjs) + 1]] <- shpMask$shp
  if (graficarObservaciones) spObjs[[length(spObjs) + 1]] <- observaciones
  if (graficarPrediccion || graficarVarianza || graficarSD) spObjs[[length(spObjs) + 1]] <- grilla
  
  xyLims <- getXYLims(spObjs=spObjs, resXImagenes, resYImagenes, ejesXYLatLong=T, factorMargen = 0)
  if (!is.null(shpMask)) { layoutSHP <- getLayoutSHP(shpMask$shp)
  } else { layoutSHP <- NULL }
  
  if (graficarObservaciones) {
    mapearPuntosGGPlot(puntos=observaciones, shpBase = shpMask$shp, escala=escalaObservaciones, 
                       nomArchResultados=changeFileExt(appendToFileName(nomArchResultados, postFijo='_Obs'), '.png'),
                       xyLims=xyLims, dibujarEjes=dibujarEjes)
  }

  if (graficarPrediccion) {
    #shpBase=shpMask$shp
    #escala=escalaObservaciones
    mapearGrillaGGPlot(grilla=grilla, shpBase=shpMask$shp, escala=escalaObservaciones,
                       nomArchResultados=changeFileExt(appendToFileName(nomArchResultados, postFijo='_Pred'), '.png'),
                       xyLims=xyLims, dibujarEscala=dibujarEscala, dibujarEjes=dibujarEjes)
  }

  if (graficarVarianza) {
    escalaVar <- escalaSD
    escalaVar$escala <- seq(0, escalaSD$escala[nrow(escalaSD)]^2, length.out=nrow(escalaSD))
    mapearGrillaGGPlot(grilla=grilla, shpBase=shpMask$shp, escala=escalaVar,
                       nomArchResultados=changeFileExt(appendToFileName(nomArchResultados, postFijo='_Var'), '.png'),
                       xyLims=xyLims, dibujarEscala=dibujarEscala, dibujarEjes=dibujarEjes, zcol=2)
  }
  
  if (graficarSD) {
    grilla$var1.sd <- sqrt(grilla@data[,2])
    mapearGrillaGGPlot(grilla=grilla, shpBase=shpMask$shp, escala=escalaVar,
                       nomArchResultados=changeFileExt(appendToFileName(nomArchResultados, postFijo='_SD'), '.png'),
                       xyLims=xyLims, dibujarEscala=dibujarEscala, dibujarEjes=dibujarEjes, zcol=3)
  }
  
  return (TRUE)
}