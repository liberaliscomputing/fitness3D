#Imports thrid-party packages for database connection, text mining, Kriging interpolation, and surface modeling
require(RODBC) 
require(tm)
require(lsa)
require(geoR) 
require(lattice)


#Extracts dataset from the local database
dataExtract <- function(method, option, duration) {
  print('Extracting data...')
  
  #Builds a query
  df.raw <- NULL
  query.select <- 'SELECT pic.patent_id, pas.org_norm_name, pas.grant_date'
  query.from <- 'FROM patent_int_classes pic, patent_assignees pas'
  if(method == 'tf-idf' || method == 'lsa') {
    query.select <- paste(c(query.select, ', pab.patent_title, pab.patent_abstract '), collapse = '')
    query.from <- paste(c(query.from, ', patent_abstracts pab '), collapse = '')
  }
  query.where <- NULL
  if(option[1] == 'ipc') {
    query.where <- paste(c('WHERE (pic.intl_class LIKE \'%', option[2], '%\') '), collapse = '')
  } else {
    for(i in 2:length(option)) {
      query.where <- paste(c(query.where, 'AND (pas.org_norm_name LIKE \'%', option[i], '%\') '), collapse = '')
    }
  }
  query.conditions <- 'AND (pic.pos = 1) AND (pic.patent_id = pas.patent_id) AND (pas.pos = 1) AND (pas.patent_id = pab.patent_id) '
  
  #Connects to the local database
  channel <- odbcConnect('MySQL_patent')
  
  #Searches rows matching with the query
  for(year in duration) {
    tmp <- sqlQuery(channel, paste(c(query.select, 
                                     query.from, 
                                     query.where, 
                                     query.conditions,
                                     'AND (pas.grant_date LIKE \'',
                                     year,
                                     '%\') ',
                                     'ORDER BY pas.grant_date ASC'), collapse = ''))
    df.raw <- rbind(df.raw, tmp)
    
  }
  rm(tmp)
  return(df.raw)
}


#Measures each patent's fitness
fitnessMeasure <- function(dataFrame, duration) {
  citation <- NULL
  #Connects to the local database
  channel <- odbcConnect('MySQL_patent')
  for(patent in dataFrame$patent_id) {
    tmp <- sqlQuery(channel, paste(c('SELECT patent_id, year, citations FROM patent_citations_by_year ',
                                     'WHERE (patent_id = ', patent, ')'), collapse = ''))
    if(length(tmp$patent_id) == 0) {  
      cnt <- 0
      citation <- c(citation, cnt)
      print(paste(c('Patent:', patent,'Fitness:', cnt), collapse = ' '))
    } else {
      cnt <- sum(subset(tmp, year <= duration[length(duration)])$citations)
      citation <- c(citation, cnt)
      print(paste(c('Patent:', patent,'Fitness:', cnt), collapse = ' '))
    }
  }
  return(citation)
}


#Pre-process data
dataPreprocess <- function(method, dataFrame) {
  if(method == 'tf-idf' || method == 'lsa') { 
    content <- within(dataFrame, merged <- paste(patent_title, patent_abstract, sep=". "))$merged
    corpus <- Corpus(VectorSource(content))
    corpus <- tm_map(corpus, tolower)
    corpus <- tm_map(corpus, removePunctuation)
    corpus <- tm_map(corpus, function(x) removeWords(x, stopwords("english")))
    corpus <- tm_map(corpus, stemDocument, language = "english")
    return(corpus)
  } else {
    bibliography <- NULL
    channel <- odbcConnect('MySQL_patent')
    for(patent in dataFrame$patent_id) {
      tmp <- NULL
      if(method == 'bc') {
        tmp <- sqlQuery(channel, paste(c('SELECT cite_to_patent_id FROM patent_citations ',
                                         'WHERE (cite_from_patent_id = ', patent, ')'), collapse = ''))
      } else if(method == 'ca') {
        tmp <- sqlQuery(channel, paste(c('SELECT cite_from_patent_id FROM patent_citations ',
                                         'WHERE (cite_to_patent_id = ', patent, ')'), collapse = ''))
      }
      citePatents <- NULL
      if(length(tmp[,1]) == 0) {
        
        citePatents <- ','
      } else {
        
        citePatents <- paste(tmp[,1], collapse = ' ')
      }
      bibliography <- c(bibliography, citePatents)
    }
    rm(tmp)
    return(bibliography)
  }
}    


#Generates a matrix
matrixGenerate <- function(preprocessedSource) {  
  matrixByMethod <- as.matrix(TermDocumentMatrix(preprocessedSource))
  return(matrixByMethod)
}


#Calculates a distance
calculateDistance <- function(method, rawMatrix) {
  if(method == 'lsa') {
    #Weighting
    mat.lsa <- lw_bintf(rawMatrix) * gw_idf(rawMatrix) 
    #Creates LSA space
    lsaSpace <- lsa(mat.lsa)
    #computes a distance matrix
    dist.mat.lsa <- dist(t(as.textmatrix(lsaSpace)))
    return(dist.mat.lsa)
  } else {
    dist.mat <- dist(t(as.matrix(td.mat)))
    return(dist.mat)
  }
}


#Projects the high dimensional space to a two-dimensional space
projectMDS <- function(distanceMatrix) {
  fit <- cmdscale(distanceMatrix, eig = TRUE, k = 2)
  points <- data.frame(x = fit$points[, 1], y = fit$points[, 2])
  return(points)
}


#Integrates data
integrateData <- function(dataFrame, coordinates, fitness) {
  df.int <- data.frame(dataFrame$patent_id, dataFrame$org_norm_name, dataFrame$grant_date, 
                       round(coordinates$x, 3), round(coordinates$y, 3), fitness)
  colnames(df.int) <- c("id", "org", "date", "x", "y", "z")
  return(df.int)
}


#Kriging interpolation
krigingIntp <- function(dataFrame) {
  #Removes duplicates in coordinates
  nonDupData <- dataFrame[!duplicated(dataFrame$x, dataFrame$y),]
  tmp <- data.frame(nonDupData$x, nonDupData$y, nonDupData$z)
  colnames(tmp) <- c("x", "y", "z")
  geoSpace <- as.geodata(tmp)
  x.range <- range(tmp[,1])
  y.range <- range(tmp[,2])
  x <- seq(from=x.range[1], to=x.range[2], by=(x.range[2] - x.range[1])/(length(tmp$x)-1))
  y <- seq(from=y.range[1], to=y.range[2], by=(y.range[2] - y.range[1])/(length(tmp$y)-1))
  xv <- rep(x,length(y))
  yv <- rep(y, each=length(x))
  in_mat <- as.matrix(cbind(xv,yv))
  q <- ksline(geoSpace, cov.model="exponential", cov.pars=c(10,3.33), nugget=0, locations=in_mat)
  intp <- data.frame(q$locations[,1], q$locations[,2], q$predict)
  colnames(intp) <- c("x", "y", "z")
  return(intp)
}


#Visualizes a surface modeled fitness landscape
visualizeSurface <- function(duration, intData, krigedData) {
  title <- NULL
  if(option[1] == 'ipc') {
    title <- paste(c('Competive Intelligence in', toupper(option[1]), toupper(option[2]), 
                     duration[1], '-', duration[length(duration)]), collapse = ' ')                
  } else {
    tmp <- NULL
    for(i in 2:length(option)) {
      tmp <- paste(c(tmp, option[i]), collapse = ' ')
    }
    title <- paste(c('Competive Intelligence in', tmp, duration[1], '-', duration[length(duration)]), collapse = ' ') 
  }
  pts <- data.frame(intData$x, intData$y, intData$z)
  colnames(pts) <- c("x", "y", "z")
  surface <- wireframe(z ~ x*y, data = krigedData, drape = T, screen = list(x = view[1], y = view[2], z = view[3]),
            xlab = "Dimension 1", ylab = "Dimension 2", zlab = "Fitness (cited times)", scales = list(arrows = F),
            main = title,
            col.regions = colorRampPalette(c("purple", "blue", "green", "yellow", "orange", "red"))(100), 
            pts = pts,
            panel.3d.wireframe =
              function(x, y, z,
                       xlim, ylim, zlim,
                       xlim.scaled, ylim.scaled, zlim.scaled,
                       pts,
                       ...) {
                panel.3dwire(x = x, y = y, z = z,
                             xlim = xlim,
                             ylim = ylim,
                             zlim = zlim,
                             xlim.scaled = xlim.scaled,
                             ylim.scaled = ylim.scaled,
                             zlim.scaled = zlim.scaled,
                             ...)
                xx <-
                  xlim.scaled[1] + diff(xlim.scaled) *
                  (pts$x - xlim[1]) / diff(xlim)
                yy <-
                  ylim.scaled[1] + diff(ylim.scaled) *
                  (pts$y - ylim[1]) / diff(ylim)
                zz <-
                  zlim.scaled[1] + diff(zlim.scaled) *
                  (pts$z - zlim[1]) / diff(zlim)
                panel.3dscatter(x = xx,
                                y = yy,
                                z = zz,
                                type = c("p", "h"), 
                                xlim = xlim,
                                ylim = ylim,
                                zlim = zlim,
                                xlim.scaled = xlim.scaled,
                                ylim.scaled = ylim.scaled,
                                zlim.scaled = zlim.scaled,
                                pch = 19, cex = 2, lty = 3,
                                col = color[as.factor(intData$org)],
                                ...)
              })
  return(surface)
}


#Defines the main fitness3D function
fitness3D <- function(method, option, color, duration, view) {
  print('Start processing...')
  start <- Sys.time()
  method <- tolower(method)
  option <- tolower(option)
  color <- tolower(color)
  
  data.extracted <- dataExtract(method, option, duration) #Data extraction
  fitness <- fitnessMeasure(data.extracted, duration) #Fitness Measurement
  preprocessed <- dataPreprocess(method, data.extracted) #Data pre-processing 
  matrixGenerated <- matrixGenerate(preprocessed) #Matrix generation
  distMatrix <- calculateDistance(method, matrixGenerated) #Distance calculation
  coordinates <- projectMDS(distMatrix) #MDS projection
  intData <- integrateData(data.extracted, coordinates, fitness) #Data integration
  kriged <- krigingIntp(intData) #Kriging interpolation
  fitnessLandscape <- visualizeSurface(duration, intData, kriged) #Fitness landscape visualization
  plot(fitnessLandscape)
  
  end <- Sys.time()
  print(paste(c('Processing times taken:', as.character(end-start)), collapse = ' '))
}