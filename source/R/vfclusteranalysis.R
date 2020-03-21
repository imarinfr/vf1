#' VF cluster analysis
#'
#' Analyse a VF for the presence of a VF defect based on the following criteria: Hodapp-Parish-Anderson 2 (HAP2), United Kington Glaucoma Treatment Study (UKGTS), Glaucoma Hemifield Test (GHT), Foster, and/or Low-pressure Glaucoma Treatment Study (LoGTS)
#' @param vf Visual fields to be analyzed as a standard vfobject
#' @param criteria select to use all criteria (default), or a single criteria: hap2, ukgts, ght, foster, or logts
#' @return results Whether VF analysis resulted in a detection of a VF defect
#' @export
vfclusteranalysis <- function( vf , criteria = "all", vf.ctrl = vfctrSunyiu24d2 )
{
  if( nrow( vf ) < 1 )
    stop("vf is empty")
  #if( length( unique( vf$tperimetry ) ) > 1 ) 
  #  stop("mixed tperimetry types: use only one perimeter test type (e.g. p24d2)")
  library(simpleboot)
  library(boot)
  library(mosaic)
  
  # obtain td, tdp, pd, pdp maps, and take out only needed columns
  #td  <- gettd(vfpwgSunyiu24d2[1,])
  td  <- gettd( vf )
  pd  <- getpd( td )

  #if( vf$tpattern[1] == "p24d2" )
  #{
    # defines indices of surrounding points to each point on vf map
    SURR_MAP <- data.frame(
      "x" =          c(-9,-3,3,9,-15,-9,-3,3,9,15,-21,-15,-9,-3,3,9,15,21,-27,-21,-15,-9,-3,3,9,15,21,-27,-21,-15,-9,-3,3,9,15,21,-21,-15,-9,-3,3,9,15,21,-15,-9,-3,3,9,15,-9,-3,3,9),
      "y" =          c(21,21,21,21,15,15,15,15,15,15,9,9,9,9,9,9,9,9,3,3,3,3,3,3,3,3,3,-3,-3,-3,-3,-3,-3,-3,-3,-3,-9,-9,-9,-9,-9,-9,-9,-9,-15,-15,-15,-15,-15,-15,-21,-21,-21,-21),
      "topleft" =    c(NA,NA,NA,NA,NA,NA,1,2,3,4,NA,NA,5,6,7,8,9,10,NA,NA,11,12,13,14,15,16,17,NA,19,20,21,22,23,24,25,26,28,29,30,31,32,33,34,35,37,38,39,40,41,42,45,46,47,48),
      "top" =        c(NA,NA,NA,NA,NA,1,2,3,4,NA,NA,5,6,7,8,9,10,NA,NA,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,29,30,31,32,33,34,35,36,38,39,40,41,42,43,46,47,48,49),
      "topright" =   c(NA,NA,NA,NA,1,2,3,4,NA,NA,5,6,7,8,9,10,NA,NA,11,12,13,14,15,16,17,18,NA,20,21,22,23,24,25,26,27,NA,30,31,32,33,34,35,36,NA,39,40,41,42,43,44,47,48,49,50),
      "left" =       c(NA,1,2,3,NA,5,6,7,8,9,NA,11,12,13,14,15,16,17,NA,19,20,21,22,23,24,25,26,NA,28,29,30,31,32,33,34,35,NA,37,38,39,40,41,42,43,NA,45,46,47,48,49,NA,51,52,53),
      "right" =      c(2,3,4,NA,6,7,8,9,10,NA,12,13,14,15,16,17,18,NA,20,21,22,23,24,25,26,27,NA,29,30,31,32,33,34,35,36,NA,38,39,40,41,42,43,44,NA,46,47,48,49,50,NA,52,53,54,NA),
      "botleft" =    c(5,6,7,8,11,12,13,14,15,16,19,20,21,22,23,24,25,26,NA,28,29,30,31,32,33,34,35,NA,NA,37,38,39,40,41,42,43,NA,NA,45,46,47,48,49,50,NA,NA,51,52,53,54,NA,NA,NA,NA),
      "bot" =        c(6,7,8,9,12,13,14,15,16,17,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,NA,37,38,39,40,41,42,43,44,NA,45,46,47,48,49,50,NA,NA,51,52,53,54,NA,NA,NA,NA,NA),
      "botright" =   c(7,8,9,10,13,14,15,16,17,18,21,22,23,24,25,26,27,NA,29,30,31,32,33,34,35,36,NA,37,38,39,40,41,42,43,44,NA,45,46,47,48,49,50,NA,NA,51,52,53,54,NA,NA,NA,NA,NA,NA),
      "refl" =       c(51,52,53,54,45,46,47,48,49,50,37,38,39,40,41,42,43,44,28,29,30,31,32,33,34,35,36,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),
      "ght.sector" = c(1,1,2,2,1,1,1,1,2,2,3,3,4,4,4,4,NA,NA,3,3,3,5,5,5,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA)
    )
    
    tdp <- gettdp( td )[which( colnames( vf ) == "l1") : which( colnames( vf ) == "l54" )]
    pdp <- getpdp( pd )[which( colnames( vf ) == "l1") : which( colnames( vf ) == "l54" )]
    td  <- td[which( colnames( vf ) == "l1") : which( colnames( vf ) == "l54" )]
    pd  <- pd[which( colnames( vf ) == "l1") : which( colnames( vf ) == "l54" )]
  #}
  #else if( vf$tpattern[1] == "p32d2" )
  #{
  #  # defines indices of surrounding points to each point on vf map
  #  SURR_MAP <- data.frame(
  #    "x" =          c(-9,-3,3,9,-15,-9,-3,3,9,15,-21,-15,-9,-3,3,9,15,21,-27,-21,-15,-9,-3,3,9,15,21,27,-27,-21,-15,-9,-3,3,9,15,21,27,-27,-21,-15,-9,-3,3,9,15,21,27,-27,-21,-15,-9,-3,3,9,15,21,27,-21,-15,-9,-3,3,9,15,21,-15,-9,-3,3,9,15,-9,-3,3,9),
  #    "y" =	         c(27,27,27,27,21,21,21,21,21,21,15,15,15,15,15,15,15,15,9,9,9,9,9,9,9,9,9,9,3,3,3,3,3,3,3,3,3,3,-3,-3,-3,-3,-3,-3,-3,-3,-3,-3,-9,-9,-9,-9,-9,-9,-9,-9,-9,-9,-15,-15,-15,-15,-15,-15,-15,-15,-21,-21,-21,-21,-21,-21,-27,-27,-27,-27),
  #    "topleft" =		 c(NA,NA,NA,NA,NA,NA,1,2,3,4,NA,NA,5,6,7,8,9,10,NA,NA,11,12,13,14,15,16,17,18,NA,19,20,21,22,23,24,25,26,27,NA,29,30,31,32,33,34,35,36,37,NA,39,40,41,42,43,44,45,46,47,49,50,51,52,53,54,55,56,59,60,61,62,63,64,67,68,69,70),
  #    "top"	=				 c(NA,NA,NA,NA,NA,1,2,3,4,NA,NA,5,6,7,8,9,10,NA,NA,11,12,13,14,15,16,17,18,NA,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,50,51,52,53,54,55,56,57,60,61,62,63,64,65,68,69,70,71),
  #    "topright"	=	 c(NA,NA,NA,NA,1,2,3,4,NA,NA,5,6,7,8,9,10,NA,NA,11,12,13,14,15,16,17,18,NA,NA,20,21,22,23,24,25,26,27,28,NA,30,31,32,33,34,35,36,37,38,NA,40,41,42,43,44,45,46,47,48,NA,51,52,53,54,55,56,57,58,61,62,63,64,65,66,69,70,71,72),
  #    "left"	=	     c(NA,1,2,3,NA,5,6,7,8,9,NA,11,12,13,14,15,16,17,NA,19,20,21,22,23,24,25,26,27,NA,29,30,31,32,33,34,35,36,37,NA,39,40,41,42,43,44,45,46,47,NA,49,50,51,52,53,54,55,56,57,NA,59,60,61,62,63,64,65,NA,67,68,69,70,71,NA,73,74,75),
  #    "right"	=      c(2,3,4,NA,6,7,8,9,10,NA,12,13,14,15,16,17,18,NA,20,21,22,23,24,25,26,27,28,NA,30,31,32,33,34,35,36,37,38,NA,40,41,42,43,44,45,46,47,48,NA,50,51,52,53,54,55,56,57,58,NA,60,61,62,63,64,65,66,NA,68,69,70,71,72,NA,74,75,76,NA),
  #    "botleft"	=    c(5,6,7,8,11,12,13,14,15,16,19,20,21,22,23,24,25,26,NA,29,30,31,32,33,34,35,36,37,NA,39,40,41,42,43,44,45,46,47,NA,49,50,51,52,53,54,55,56,57,NA,NA,59,60,61,62,63,64,65,66,NA,NA,67,68,69,70,71,72,NA,NA,73,74,75,76,NA,NA,NA,NA),
  #    "bot"	=        c(6,7,8,9,12,13,14,15,16,17,20,21,22,23,24,25,26,27,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,NA,59,60,61,62,63,64,65,66,NA,NA,67,68,69,70,71,72,NA,NA,73,74,75,76,NA,NA,NA,NA,NA),			
  #    "botright" =   c(7,8,9,10,13,14,15,16,17,18,21,22,23,24,25,26,27,28,30,31,32,33,34,35,36,37,38,NA,40,41,42,43,44,45,46,47,48,NA,50,51,52,53,54,55,56,57,58,NA,59,60,61,62,63,64,65,66,NA,NA,67,68,69,70,71,72,NA,NA,73,74,75,76,NA,NA,NA,NA,NA,NA),
  #    "refl" =       c(73,74,75,76,67,68,69,70,71,72,59,60,61,62,63,64,65,66,49,50,51,52,53,54,55,56,57,58,39,40,41,42,43,44,45,46,47,48,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),
  #    "ght.sector" = c(NA,NA,NA,NA,NA,1,1,2,2,NA,NA,1,1,1,1,2,2,NA,NA,3,3,4,4,4,4,NA,NA,NA,3,3,3,5,5,5,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA)
  #    
  #  )		
  #  
  #  tdp <- tdpval( td )[which( colnames( vf ) == "l1" ) : which( colnames( vf ) == "l76" )]
  #  pdp <- getpdp( pd )[which( colnames( vf ) == "l1" ) : which( colnames( vf ) == "l76" )]
  #  td  <- td[which( colnames( vf ) == "l1" ) : which( colnames( vf ) == "l76" )]
  #  pd  <- pd[which( colnames( vf ) == "l1" ) : which( colnames( vf ) == "l76" )]
  #}
  #else
  #{
  #  stop( "invalid tperimtery test type: pass only a p24d2 or p32d2 test type" )
  #}
 
  if( ( criteria == "all") || ( criteria == "ght") )
  { 
    #' Calculate limits of normality as described in Asman & Heijl, Arch Opthalmol, 1992
    #' 1. get about 200+ vfs of normal eyes, preferably 1 unique eye per patient
    #' 2. select 500 random samples, each with 15 vfs, from the total 200+ vfs, making sure each sample has 1 unique eye per patient
    #' 3. use the bootstrap method to find confidence limits for gh, sum.sup, sum.inf, up.down using the 500 random samples

    CTRL_VF_NUM = nrow(vf.ctrl)
    SUMS_NUM = 3
    GHT_SECTORS_NUM = 5
    SECTORS_LIM_NUM = 6
    GH_LIM_NUM = 2
    SAMPLE_SIZE = 25
    SAMPLE_NUM = 500
    
    td.ctrl  <- gettd( vf.ctrl )
    pd.ctrl  <- getpd( td.ctrl )
    pdp.ctrl <- getpdp( pd.ctrl )[which( colnames( vf.ctrl ) == "l1") : which( colnames( vf.ctrl ) == "l54" )]
    pd.ctrl  <- pd.ctrl[which( colnames( vf.ctrl ) == "l1") : which( colnames( vf.ctrl ) == "l54" )]
    
    # data structures
    sector.sums <- array(NA, c( GHT_SECTORS_NUM, SUMS_NUM, CTRL_VF_NUM ) )
    colnames( sector.sums ) <- c( "sum.sup", "sum.inf", "up.down" )
    
    gh <- array(NA, c( 1, 1, CTRL_VF_NUM ) )
    colnames( gh ) <- c( "gh" )
    
    for( v in 1:CTRL_VF_NUM )
    {
      scores <- pdp.ctrl[v,]
      for( i in 1:length( pdp.ctrl[v,] ) )
      {
        if( !( is.na( pdp.ctrl[v,i] ) ) )
        {
          if( pdp.ctrl[v,i] > 5 )
            scores[i] <- 0
          else if( pdp.ctrl[v,i] > 2 )
            scores[i] <- 2
          else if( pdp.ctrl[v,i] > 1 )
            scores[i] <- 5
          else
            scores[i] <- 10 * abs( pd.ctrl[v,i] )      
        }
      }
      for( sector in 1:GHT_SECTORS_NUM )
      {
        sector.sums[ sector,"sum.sup",v ] <- sum( scores[ which( SURR_MAP$ght.sector == sector ) ] )
        sector.sums[ sector,"sum.inf",v ] <- sum( scores[ SURR_MAP[ which( SURR_MAP$ght.sector == sector ), "refl" ] ] )
        sector.sums[ sector,"up.down",v ] <- sector.sums[ sector,"sum.sup",v ] - sector.sums[ sector,"sum.inf",v ]
      }
      
      gh[ 1,"gh", ] <- getgh( td.ctrl )
    }
    
    sector.lims <- data.frame( array( NA, c( GHT_SECTORS_NUM, SECTORS_LIM_NUM ) ) )
    colnames( sector.lims ) <- c("sum.sup99.5", "sum.inf99.5", "up.down0.5", "up.down99.5", "up.down1.5", "up.down98.5" )
    
    gh.lims <- data.frame( array( NA, c( 1, GH_LIM_NUM ) ) )
    colnames( gh.lims ) <- c( "gh0.5", "gh99.5" )
    
    for( sector in 1:GHT_SECTORS_NUM )
    {
      #sector.lims[ sector,"sum.sup99.5"] <- as.double( boot.ci( one.boot( sector.sums[ sector,"sum.sup", ], mean, R = 500 ), conf = 0.995, type = "norm")[[4]][3] )
      #sector.lims[ sector,"sum.inf99.5"] <- as.double( boot.ci( one.boot( sector.sums[ sector,"sum.inf", ], mean, R = 500 ), conf = 0.995, type = "norm")[[4]][2] )
      #sector.lims[ sector,3:4] <- as.double( boot.ci( one.boot( sector.sums[ sector,"up.down", ], mean, R = 500 ), conf = 0.995, type = "norm")[[4]][2:3] )
      #sector.lims[ sector,5:6] <- as.double( boot.ci( one.boot( sector.sums[ sector,"up.down", ], mean, R = 500 ), conf = 0.985, type = "norm")[[4]][2:3] )    
      
      # generate bootstrapped confidence limits for each sector sum and up-down
      bs.lims.sum.sup   <- do(SAMPLE_NUM) * suppressMessages( confint( resample( sector.sums[ sector,"sum.sup", ], SAMPLE_SIZE), level = 0.99, method = "quantile" ) )
      bs.lims.sum.inf   <- do(SAMPLE_NUM) * suppressMessages( confint( resample( sector.sums[ sector,"sum.inf", ], SAMPLE_SIZE), level = 0.99, method = "quantile" ) )
      bs.lims.up.down99 <- do(SAMPLE_NUM) * suppressMessages( confint( resample( sector.sums[ sector,"up.down", ], SAMPLE_SIZE), level = 0.99, method = "quantile" ) )
      bs.lims.up.down97 <- do(SAMPLE_NUM) * suppressMessages( confint( resample( sector.sums[ sector,"up.down", ], SAMPLE_SIZE), level = 0.97, method = "quantile" ) )
      
      sector.lims[ sector,"sum.sup99.5"] <- mean( bs.lims.sum.sup$X99.5. )
      sector.lims[ sector,"sum.inf99.5"] <- mean( bs.lims.sum.inf$X99.5. )                
      sector.lims[ sector,"up.down0.5"]  <- mean( bs.lims.up.down99$X0.5. )
      sector.lims[ sector,"up.down99.5"] <- mean( bs.lims.up.down99$X99.5. )
      sector.lims[ sector,"up.down1.5"]  <- mean( bs.lims.up.down97$X1.5. )
      sector.lims[ sector,"up.down98.5"] <- mean( bs.lims.up.down97$X98.5. ) 
      }
    
    #gh.lims[1,1:2] <- as.double( data.frame( boot.ci( one.boot( gh[ 1,"gh", ], mean, R = 500 ), conf = 0.995, type = "norm")[4] )[2:3] )
    bs.lims.gh <- do(SAMPLE_NUM) * suppressMessages( confint( resample( gh[ 1,"gh", ], SAMPLE_SIZE), level = 0.99, method = "quantile" ) )
    gh.lims$gh99.5 <- mean( bs.lims.gh$X99.5. )
    gh.lims$gh0.5 <- mean( bs.lims.gh$X0.5. )
    
    ght.lims = list( "gh" = gh.lims, "sector" = sector.lims )
    
    assign("sector.sums", sector.sums, envir = globalenv())
    assign("cgh", gh, envir = globalenv())
    assign("ght.lims", ght.lims, envir = globalenv())
    #print(ght.lims)
  }

  # implement results data frame  
  if( criteria == "all")
  {
    result <- data.frame( array( FALSE, c( nrow( vf ), 5 ) ) )
    colnames( result ) = c( "hap2", "ukgts", "ght", "foster", "logts" )
  }
  else
  {
    result <- data.frame( array( FALSE, c( nrow( vf ), 1 ) ) )
    colnames( result ) <- criteria
  }
  
  rownames( result ) <- rownames( vf )
  
  # analyze vf for defect
  for( i in 1:nrow(vf_data) )
  {
    if( criteria == "all" )
      result[i,] <- c( hap2( vf[i,], td[i,], pd[i,], pdp[i,], ght.lims, SURR_MAP ), 
                       ukgts( vf[i,], td[i,], tdp[i,], SURR_MAP ), 
                       ght( vf[i,], td[i,], pd[i,], pdp[i,], ght.lims, SURR_MAP ), 
                       foster( vf[i,], td[i,], pd[i,], pdp[i,], ght.lims, SURR_MAP ), 
                       logts( td[i,], SURR_MAP ) )
    
    else if( criteria == "hap2" )
      result[i,] <- hap2( vf[i,], td[i,], pd[i,], pdp[i,], ght.lims, SURR_MAP )
    
    else if( criteria == "ukgts" )
      result[i,] <- ukgts( vf[i,], td[i,], tdp[i,], SURR_MAP )
    
    else if( criteria == "ght" )
    {
      result[i,] <- ght( vf[i,], td[i,], pd[i,], pdp[i,], ght.lims, SURR_MAP )
    }
    else if( criteria == "foster" )
      result[i,] <- foster( vf[i,], td[i,], pd[i,], pdp[i,], ght.lims, SURR_MAP )
    
    else if( criteria == "logts" )
      result[i,] <- logts( td[i,], SURR_MAP )
    
    else
      stop( "invalid criteria: select all (default), hap2, ukgts, ght, foster, or logts" )
  }
  
  return( result )
}

vfclusterverification <- function(vf = vf_d, SAMPLE_SIZE = 10000, criteria = "all", vf.ctrl = vfctrSunyiu24d2 )
{
  library(simpleboot)
  library(boot)
  library(mosaic)
  
  #vf_data = read.csv("GitHub/vf1/source/R/vf_data.csv")
  #crit_res = read.csv("GitHub/vf1/source/R/crit1.csv")
  ## obtain td, tdp, pd, pdp maps, and take out only needed columns
  #vf_d  <- vf_data[which( colnames( vf_data ) == "vf.1")  : which( colnames( vf_data ) == "vf.54" )]
  #dummy_vf = vfctrSunyiu24d2[1,]
  #for( i in 1:(nrow(vf_data)-1))
  #{
  #  dummy_vf = rbind(dummy_vf, vfctrSunyiu24d2[1,])
  #}
  #dummy_vf[1:nrow(dummy_vf), 11:64] = vf_d
  #vf_d = dummy_vf
  #rm(dummy_vf)
  #vf_d$eye = "OD"
  
  td  <- vf_data[which( colnames( vf_data ) == "td.1")  : which( colnames( vf_data ) == "td.54" )]
  tdp <- vf_data[which( colnames( vf_data ) == "tdp.1") : which( colnames( vf_data ) == "tdp.54" )]
  pd  <- vf_data[which( colnames( vf_data ) == "pd.1")  : which( colnames( vf_data ) == "pd.54" )]
  pdp <- vf_data[which( colnames( vf_data ) == "pdp.1") : which( colnames( vf_data ) == "pdp.54" )]
  
  colnames(td) = colnames(vf)[11:64]
  colnames(pd) = colnames(vf)[11:64]
  colnames(tdp) = colnames(vf)[11:64]
  colnames(pdp) = colnames(vf)[11:64]
  
  
  # defines indices of surrounding points to each point on vf map
  SURR_MAP <- data.frame(
    "x" =          c(-9,-3,3,9,-15,-9,-3,3,9,15,-21,-15,-9,-3,3,9,15,21,-27,-21,-15,-9,-3,3,9,15,21,-27,-21,-15,-9,-3,3,9,15,21,-21,-15,-9,-3,3,9,15,21,-15,-9,-3,3,9,15,-9,-3,3,9),
    "y" =          c(21,21,21,21,15,15,15,15,15,15,9,9,9,9,9,9,9,9,3,3,3,3,3,3,3,3,3,-3,-3,-3,-3,-3,-3,-3,-3,-3,-9,-9,-9,-9,-9,-9,-9,-9,-15,-15,-15,-15,-15,-15,-21,-21,-21,-21),
    "topleft" =    c(NA,NA,NA,NA,NA,NA,1,2,3,4,NA,NA,5,6,7,8,9,10,NA,NA,11,12,13,14,15,16,17,NA,19,20,21,22,23,24,25,26,28,29,30,31,32,33,34,35,37,38,39,40,41,42,45,46,47,48),
    "top" =        c(NA,NA,NA,NA,NA,1,2,3,4,NA,NA,5,6,7,8,9,10,NA,NA,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,29,30,31,32,33,34,35,36,38,39,40,41,42,43,46,47,48,49),
    "topright" =   c(NA,NA,NA,NA,1,2,3,4,NA,NA,5,6,7,8,9,10,NA,NA,11,12,13,14,15,16,17,18,NA,20,21,22,23,24,25,26,27,NA,30,31,32,33,34,35,36,NA,39,40,41,42,43,44,47,48,49,50),
    "left" =       c(NA,1,2,3,NA,5,6,7,8,9,NA,11,12,13,14,15,16,17,NA,19,20,21,22,23,24,25,26,NA,28,29,30,31,32,33,34,35,NA,37,38,39,40,41,42,43,NA,45,46,47,48,49,NA,51,52,53),
    "right" =      c(2,3,4,NA,6,7,8,9,10,NA,12,13,14,15,16,17,18,NA,20,21,22,23,24,25,26,27,NA,29,30,31,32,33,34,35,36,NA,38,39,40,41,42,43,44,NA,46,47,48,49,50,NA,52,53,54,NA),
    "botleft" =    c(5,6,7,8,11,12,13,14,15,16,19,20,21,22,23,24,25,26,NA,28,29,30,31,32,33,34,35,NA,NA,37,38,39,40,41,42,43,NA,NA,45,46,47,48,49,50,NA,NA,51,52,53,54,NA,NA,NA,NA),
    "bot" =        c(6,7,8,9,12,13,14,15,16,17,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,NA,37,38,39,40,41,42,43,44,NA,45,46,47,48,49,50,NA,NA,51,52,53,54,NA,NA,NA,NA,NA),
    "botright" =   c(7,8,9,10,13,14,15,16,17,18,21,22,23,24,25,26,27,NA,29,30,31,32,33,34,35,36,NA,37,38,39,40,41,42,43,44,NA,45,46,47,48,49,50,NA,NA,51,52,53,54,NA,NA,NA,NA,NA,NA),
    "refl" =       c(51,52,53,54,45,46,47,48,49,50,37,38,39,40,41,42,43,44,28,29,30,31,32,33,34,35,36,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),
    "ght.sector" = c(1,1,2,2,1,1,1,1,2,2,3,3,4,4,4,4,NA,NA,3,3,3,5,5,5,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA)
  )
  
  if( ( criteria == "all") || ( criteria == "ght") )
  { 
    #' Calculate limits of normality as described in Asman & Heijl, Arch Opthalmol, 1992
    #' 1. get about 200+ vfs of normal eyes, preferably 1 unique eye per patient
    #' 2. select 500 random samples, each with 15 vfs, from the total 200+ vfs, making sure each sample has 1 unique eye per patient
    #' 3. use the bootstrap method to find confidence limits for gh, sum.sup, sum.inf, up.down using the 500 random samples
    
    CTRL_VF_NUM = nrow(vf.ctrl)
    SUMS_NUM = 3
    GHT_SECTORS_NUM = 5
    SECTORS_LIM_NUM = 6
    GH_LIM_NUM = 2
    #SAMPLE_SIZE = 10000
    SAMPLE_NUM = 500
    
    td.ctrl  <- gettd( vf.ctrl )
    pd.ctrl  <- getpd( td.ctrl )
    pdp.ctrl <- getpdp( pd.ctrl )[which( colnames( vf.ctrl ) == "l1") : which( colnames( vf.ctrl ) == "l54" )]
    pd.ctrl  <- pd.ctrl[which( colnames( vf.ctrl ) == "l1") : which( colnames( vf.ctrl ) == "l54" )]
    
    # data structures
    sector.sums <- array(NA, c( GHT_SECTORS_NUM, SUMS_NUM, CTRL_VF_NUM ) )
    colnames( sector.sums ) <- c( "sum.sup", "sum.inf", "up.down" )
    
    gh <- array(NA, c( 1, 1, CTRL_VF_NUM ) )
    colnames( gh ) <- c( "gh" )
    
    for( v in 1:CTRL_VF_NUM )
    {
      scores <- pdp.ctrl[v,]
      for( i in 1:length( pdp.ctrl[v,] ) )
      {
        if( !( is.na( pdp.ctrl[v,i] ) ) )
        {
          if( pdp.ctrl[v,i] > 5 )
            scores[i] <- 0
          else if( pdp.ctrl[v,i] > 2 )
            scores[i] <- 2
          else if( pdp.ctrl[v,i] > 1 )
            scores[i] <- 5
          else
            scores[i] <- 10 * abs( pd.ctrl[v,i] / sunyiu_24d2$luts$pd["1%",i] )   
        }
      }
      for( sector in 1:GHT_SECTORS_NUM )
      {
        sector.sums[ sector,"sum.sup",v ] <- sum( scores[ which( SURR_MAP$ght.sector == sector ) ] )
        sector.sums[ sector,"sum.inf",v ] <- sum( scores[ SURR_MAP[ which( SURR_MAP$ght.sector == sector ), "refl" ] ] )
        sector.sums[ sector,"up.down",v ] <- sector.sums[ sector,"sum.sup",v ] - sector.sums[ sector,"sum.inf",v ]
      }
      
      gh[ 1,"gh", ] <- getgh( td.ctrl )
    }
    sector.lims <- data.frame( array( NA, c( GHT_SECTORS_NUM, SECTORS_LIM_NUM ) ) )
    colnames( sector.lims ) <- c("sum.sup99.5", "sum.inf99.5", "up.down0.5", "up.down99.5", "up.down1.5", "up.down98.5" )
    
    gh.lims <- data.frame( array( NA, c( 1, GH_LIM_NUM ) ) )
    colnames( gh.lims ) <- c( "gh0.5", "gh99.5" )
    
    for( sector in 1:GHT_SECTORS_NUM )
    {
      #sector.lims[ sector,"sum.sup99.5"] <- as.double( boot.ci( one.boot( sector.sums[ sector,"sum.sup", ], mean, R = SAMPLE_NUM ), conf = 0.995, type = "norm")[[4]][3] )
      #sector.lims[ sector,"sum.inf99.5"] <- as.double( boot.ci( one.boot( sector.sums[ sector,"sum.inf", ], mean, R = SAMPLE_NUM ), conf = 0.995, type = "norm")[[4]][2] )
      #sector.lims[ sector,3:4] <- as.double( boot.ci( one.boot( sector.sums[ sector,"up.down", ], mean, R = SAMPLE_NUM ), conf = 0.995, type = "norm")[[4]][2:3] )
      #sector.lims[ sector,5:6] <- as.double( boot.ci( one.boot( sector.sums[ sector,"up.down", ], mean, R = SAMPLE_NUM ), conf = 0.985, type = "norm")[[4]][2:3] )    
      
      # generate bootstrapped confidence limits for each sector sum and up-down difference
      bs.lims.sum.sup   <- do(SAMPLE_NUM) * suppressMessages( confint( resample( sector.sums[ sector,"sum.sup", ], SAMPLE_SIZE), level = 0.99, method = "quantile" ) )
      bs.lims.sum.inf   <- do(SAMPLE_NUM) * suppressMessages( confint( resample( sector.sums[ sector,"sum.inf", ], SAMPLE_SIZE), level = 0.99, method = "quantile" ) )
      bs.lims.up.down99 <- do(SAMPLE_NUM) * suppressMessages( confint( resample( sector.sums[ sector,"up.down", ], SAMPLE_SIZE), level = 0.98, method = "quantile" ) )
      bs.lims.up.down97 <- do(SAMPLE_NUM) * suppressMessages( confint( resample( sector.sums[ sector,"up.down", ], SAMPLE_SIZE), level = 0.94, method = "quantile" ) )
      
      sector.lims[ sector,"sum.sup99.5"] <- mean( bs.lims.sum.sup$X99.5. )
      sector.lims[ sector,"sum.inf99.5"] <- mean( bs.lims.sum.inf$X99.5. )                
      sector.lims[ sector,"up.down0.5"]  <- mean( bs.lims.up.down99$X1. )
      sector.lims[ sector,"up.down99.5"] <- mean( bs.lims.up.down99$X99. )
      sector.lims[ sector,"up.down1.5"]  <- mean( bs.lims.up.down97$X3. )
      sector.lims[ sector,"up.down98.5"] <- mean( bs.lims.up.down97$X97. ) 
    }
    #gh.lims[1,1:2] <- boot.ci( one.boot( gh[ 1,"gh", ], mean, R = SAMPLE_NUM ), conf = 0.995, type = "norm")[[4]][2:3]
    bs.lims.gh <- do(SAMPLE_NUM) * suppressMessages( confint( resample( gh[ 1,"gh", ], SAMPLE_SIZE), level = 0.99, method = "quantile" ) )
    gh.lims$gh99.5 <- mean( bs.lims.gh$X99.5. )
    gh.lims$gh0.5 <- mean( bs.lims.gh$X0.5. )
    
    ght.lims = list( "gh" = gh.lims, "sector" = sector.lims )
    
    assign("ss", sector.sums, envir = globalenv())
    assign("g", gh, envir = globalenv())
    assign("ght.lims", ght.lims, envir = globalenv())
    print(ght.lims)
  }
  
  # implement results data frame
  if( criteria == "all")
  {
    result <- data.frame( array( FALSE, c( nrow( vf ), 5 ) ) )
    colnames( result ) = c( "hap2", "ukgts", "ght", "foster", "logts" )
  }
  else
  {
    result <- data.frame( array( FALSE, c( nrow( vf ), 1 ) ) )
    colnames( result ) <- criteria
  }
  
  rownames( result ) <- rownames( vf )
  
  # analyze vf for defect
  print("ok")
  for( i in 1:nrow(vf) )
  {
    #print(i)
    if( criteria == "all" )
      result[i,] <- c( hap2( vf[i,], td[i,], pd[i,], pdp[i,], ght.lims, SURR_MAP ), 
                       ukgts( vf[i,], td[i,], tdp[i,], SURR_MAP ), 
                       ght( vf[i,], td[i,], pd[i,], pdp[i,], ght.lims, SURR_MAP ), 
                       foster( vf[i,], td[i,], pd[i,], pdp[i,], ght.lims, SURR_MAP ), 
                       logts( td[i,], SURR_MAP ) )
    
    else if( criteria == "hap2" )
      result[i,] <- hap2( vf[i,], td[i,], pd[i,], pdp[i,], ght.lims, SURR_MAP )
    
    else if( criteria == "ukgts" )
      result[i,] <- ukgts( vf[i,], td[i,], tdp[i,], SURR_MAP )
    
    else if( criteria == "ght" )
    {
      result[i,] <- ght( vf[i,], td[i,], pd[i,], pdp[i,], ght.lims, SURR_MAP )
      #if( ght( vf[i,], pd[i,], pdp[i,], ght.lims, SURR_MAP ) == "Outside normal limits")
      #  result[i,] <- TRUE
      #else
      #  result[i,] <- FALSE
    }
    else if( criteria == "foster" )
      result[i,] <- foster( vf[i,], td[i,], pd[i,], pdp[i,], ght.lims, SURR_MAP )
    
    else if( criteria == "logts" )
      result[i,] <- logts( td[i,], SURR_MAP )
    
    else
      stop( "invalid criteria: select all (default), hap2, ukgts, ght, foster, or logts" )
  }
  
  analysis = data.frame("n" = c(length(which(result$ght=="Outside normal limits")),
                                length(which(result$ght=="Within normal limits")),
                                length(which(result$ght=="Abnormally high sensitivity")),
                                length(which(result$ght=="Borderline")),
                                length(which(result$ght=="Borderline and general reduction in sensitivity")),
                                length(which(result$ght=="General reduction in sensitivity"))
                                ),
                        "p" = c(length(which(vf_data[which(vf_data$ght=="Outside normal limits"),"ght"] == result[which(vf_data$ght=="Outside normal limits"),"ght"]))/length(which(vf_data$ght=="Outside normal limits")),
                                length(which(vf_data[which(vf_data$ght=="Within normal limits"),"ght"] == result[which(vf_data$ght=="Within normal limits"),"ght"]))/length(which(vf_data$ght=="Within normal limits")),
                                length(which(vf_data[which(vf_data$ght=="Abnormally high sensitivity"),"ght"] == result[which(vf_data$ght=="Abnormally high sensitivity"),"ght"]))/length(which(vf_data$ght=="Abnormally high sensitivity")),
                                length(which(vf_data[which(vf_data$ght=="Borderline"),"ght"] == result[which(vf_data$ght=="Borderline"),"ght"]))/length(which(vf_data$ght=="Borderline")),
                                length(which(vf_data[which(vf_data$ght=="Borderline and general reduction in sensitivity"),"ght"] == result[which(vf_data$ght=="Borderline and general reduction in sensitivity"),"ght"]))/length(which(vf_data$ght=="Borderline and general reduction in sensitivity")),
                                length(which(vf_data[which(vf_data$ght=="General reduction in sensitivity"),"ght"] == result[which(vf_data$ght=="General reduction in sensitivity"),"ght"]))/length(which(vf_data$ght=="General reduction in sensitivity"))
                                )
                        )
  print(analysis)
  
  return( result )
}

#' Hoddap-Parrish-Anderson 2 criteria (HAP2)
#'
#' Analyse a VF for the presence of a VF defect based on the following criteria: GHT outside normal limits; OR cluster of 3 points on the pattern deviation plot depressed at p<5%, one of which depressed at p<1%; OR PSD with p<5%
#' @param vf visual fields to be analyzed as a standard vfovject
#' @param pd pattern deviation map to be analyzed
#' @param pdp pattern deviation probability map to be analyzed
#' @param lims limits of normality for the GHT criteria
#' @param helper helper map tht specifies the surrounding points of each point
#' @return BOOL Whether VF analysis resulted in a detection of a visual field defect
hap2 = function( vf, pd, pdp, lims, helper )
{
  print( "hap2")
  # check PSD result of vf
  if( getglp( getgl( vf ) )["psd"] <= 5 )
    return( TRUE )  
  
  # check GHT result of vf
  if( "Outside normal limits" == ght( vf, pd, pdp, lims, helper ) )
    return( FALSE )
  
  # get index of first and last surround map points
  m_i <- which( colnames( helper ) == "topleft" )
  m_f <- which( colnames( helper ) == "botright" )
  
  # check cluster of 3 points
  # start by checking each point i in the pdp
  for( i in 1:length( pdp ) )
  {
    if( !( is.na( pdp[i] ) ) )
    {
      if( pdp[i] <= 5 )
      {
        print(paste0("i=", i))
        
        # then, check all surrounding points for 2nd point j of contiguous defect
        for( m in m_i:m_f )
        {
          j <- helper[i,m]
          
          if( !( is.na( j ) ) && !( is.na( pdp[j] ) ) )
          {
            if( pdp[j] <= 5 )
            {
              print(paste0("j=", j))
              
              # then, check all surrounding points for 3rd point k of contiguous defect
              for( n in m_i:m_f )
              {
                k <- helper[j,n]
                
                if( !( is.na( k ) ) && !( is.na( pdp[k] ) ) )
                {
                  if( pdp[k] <= 5 )
                  {
                    # check that 1st and 3rd points are not the same
                    if( i != k )
                    {
                      # check if at least one point is less than 1%
                      if( ( pdp[i] <= 1 ) || ( pdp[j] <= 1 ) || ( pdp[k] <= 1 ) )
                      {
                        print(paste0("k=", k))
                        
                        return( TRUE )
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  return( FALSE )
}

#' United Kingdom Glaucoma Treatment Study criteria (UKGTS)
#'
#' Analyse a VF for the presence of a VF defect based on the following criteria: 2 or more contiguous points with P<0.01 loss or more; OR 3 or more contiguous points with P<0.05 loss or more; OR a 10-dB difference across the nasal horizontal midline at 2 or more adjacent points in the total deviation plot.
#' @param vf visual field map to be analyzed
#' @param td total deviation map to be analyzed
#' @param tdp total deviation probability map to be analyzed
#' @param helper helper map tht specifies the surrounding points of each point
#' @return BOOL Whether VF analysis resulted in a detection of a visual field defect
ukgts <- function( vf, td, tdp, helper )
{
  print( "ukgts")
  
  # get index of first and last surround map points
  m_i <- which( colnames( helper ) == "topleft" )
  m_f <- which( colnames( helper ) == "botright" )
  
  # check cluster of 2 points
  # start by checking each point i in the tdp
  for( i in 1:length( tdp ) )
  {
    if( !( is.na( tdp[i] ) ) )
    {
      if( tdp[i] <= 1 )
      {
        print(paste0("i=", i))
        
        # then, check all surrounding points for 2nd point j of contiguous defect
        for( m in m_i:m_f )
        {
          j <- helper[i,m]
          
          if( !( is.na( j ) ) && !( is.na( tdp[j] ) ) )
          {
            if( tdp[j] <= 1 )
            {
              print(paste0("j=", j))
              
              return( TRUE )
            }
          }
        }
      }
    }
  }
  
  # check cluster of 3 points
  # start by checking each point i in the tdp
  for( i in 1:length( tdp ) )
  {
    if( !( is.na( tdp[i] ) ) )
    {
      if( tdp[i] <= 5 )
      {
        print(paste0("i=", i))
        
        # then, check all surrounding points for 2nd point j of contiguous defect
        for( m in m_i:m_f )
        {
          j <- helper[i,m]
          
          if( !( is.na( j ) ) && !( is.na( tdp[j] ) ) )
          {
            if( tdp[j] <= 5 )
            {
              print(paste0("j=", j))
              
              # then, check all surrounding points for 3rd point k of contiguous defect
              for( n in m_i:m_f )
              {
                k <- helper[j,n]
                
                if( !( is.na( k ) ) && !( is.na( tdp[k] ) ) )
                {
                  if( tdp[k] <= 5 )
                  {
                    # check that 1st and 3rd points are not the same
                    if( i != k )
                    {
                      print(paste0("k=", k))
                      
                      return( TRUE )
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  
  # get index of first and last reflection points
  #if( vf$tpattern == "p24d2" )
  #{
    ir_i <- which( colnames( td ) == "l19" )
    ir_f <- which( colnames( td ) == "l23" )
  #}
  #else
  #{
  #  ir_i <- which( colnames( td ) == "l29" )
  #  ir_f <- which( colnames( td ) == "l33" )    
  #}

  
  # check reflection points
  # start by checking each point i in the td
  for( i in ir_i:ir_f )
  {
    ir <- helper[i,"refl"]
    
    if( !( is.na( td[i] ) ) && !( is.na( td[ir] ) ) )
    {
      if( abs( td[i] - td[ir] ) >= 5 )
      {
        print(paste0("i=", i, " ir=", ir))
        
        for( m in m_i:m_f )
        {
          j <- helper[i,m]
          
          jr <- helper[j,"refl"]
          
          if( !( is.na( j ) ) && !( is.na( jr ) ) )
          {
            if( abs( td[j] - td[jr] ) >= 5 )
            {
              if( ( i != j ) && ( ir != j ) )
              {
                print(paste0("j=", j, " jr=", jr))
                
                return( TRUE )
              }
            }
          }
        }
      }
    }
  }
  return(FALSE)
}

#' Glaucoma Hemifield Test criteria (GHT)
#'
#' Analyse a VF for the presence of a VF defect based on the following criteria: Glaucoma Hemifield Test (GHT) see Asman & Heijl Arch Ophthalmol - vol 110 1992
#' @param vf visual fields to be analyzed as a standard vfovject
#' @param pd pattern deviation map to be analyzed
#' @param pdp pattern deviation probability map to be analyzed
#' @param lims limits of normality required for the GHT algorithm
#' @param helper helper map tht specifies the surrounding points of each point
#' @return BOOL Whether VF analysis resulted in a detection of a visual field defect
ght <- function( vf, td, pd, pdp, lims, helper )
{
  # calculation of general height gh of visual field
  gh <- getgh( cbind( vf[1,1:10], td ) )

  if( gh > lims$gh$gh99.5 )
  {
    return( "Abnormally high sensitivity" )
  }
  # calculation of point-by-point probability scores

  score <- pdp
  for( i in 1:length( pdp ) )
  {
    if( is.na( pdp[i] ) )
      score[i] <- 0
    else if( pdp[i] > 5)
      score[i] <- 0
    else if( pdp[i] > 2)
      score[i] <- 2
    else if( pdp[i] > 1)
      score[i] <- 5
    else
      score[i] <- 10 * abs( pd[i] / sunyiu_24d2$luts$pd["1%",i] )      
  }
  #  print(score)
  
  # calculation of the 10 sums of scores from the 10 sectors, and the 5 up-down sector differences
  sums <- data.frame( matrix( NA, nrow = 5, ncol = 3 ))
  colnames( sums ) <- c( "sum.sup", "sum.inf", "up.down")
  for( sector in 1:5 )
  {
    sums[sector,"sum.sup"] <- sum( score[ which( helper$ght.sector == sector ) ] )
    sums[sector,"sum.inf"] <- sum( score[ helper[ which( helper$ght.sector == sector ), "refl" ] ] )
    sums[sector,"up.down"] <- sums[sector,"sum.sup"] - sums[sector,"sum.inf"]
  }
  #print(gh)
  #print(sums)
  if( ( length( which( ( sums$up.down > lims$sector$up.down99.5 ) == TRUE ) ) > 0 ) ||
      ( length( which( ( sums$up.down < lims$sector$up.down0.5 ) == TRUE ) ) > 0 ) )
  {
    return( "Outside normal limits" )
  }

  if( ( length( which( ( sums$sum.sup > lims$sector$sum.sup99.5 ) == TRUE ) ) > 0 ) ||
      ( length( which( ( sums$sum.inf > lims$sector$sum.inf99.5 ) == TRUE ) ) > 0 ) )
  {
    return( "Outside normal limits" )
  }

  if( ( length( which( ( sums$up.down > lims$sector$up.down98.5 ) == TRUE ) ) > 0 ) ||
      ( length( which( ( sums$up.down < lims$sector$up.down1.5 ) == TRUE ) ) > 0 ) )
  {
    if( gh < lims$gh$gh0.5 )
      return( "Borderline and general reduction in sensitivity" )
    else
      return( "Borderline" )
  }

  if( gh < lims$gh$gh0.5 )
    return( "General reduction in sensitivity" )
  
  return( "Within normal limits")
}

#' Foster criteria (FOST)
#'
#' Analyse a VF for the presence of a VF defect based on the following criteria: Glaucoma Hemifield Test (GHT) "outside normal limits"; AND a cluster of three contiguous points at the 5% level on the pattern deviation plot
#' @param vf visual fields to be analyzed as a standard vfovject
#' @param pd pattern deviation map to be analyzed
#' @param pdp pattern deviation probability map to be analyzed
#' @param helper helper map tht specifies the surrounding points of each point
#' @return BOOL Whether VF analysis resulted in a detection of a visual field defect
foster <- function( vf, pd, pdp, lims, helper )
{
  print( "foster" )
  
  # obtain GHT result of vf
  if( "Outside normal limits" == ght( vf, pd, pdp, lims, helper ) )
    return( FALSE )
  
  # get index of first and last surround map points
  m_i <- which( colnames( helper ) == "topleft" )
  m_f <- which( colnames( helper ) == "botright" )
  
  # start by checking each point i in the pdp
  for( i in 1:length( pdp ) )
  {
    if( !( is.na( pdp[i] ) ) )
    {
      if( pdp[i] <= 5 )
      {
        print(paste0("i=", i))
        
        # then, check all surrounding points for 2nd point j of contiguous defect
        for( m in m_i:m_f )
        {
          j <- helper[i,m]
          
          if( !( is.na( j ) ) && !( is.na( pdp[j] ) ) )
          {
            if( pdp[j] <= 5 )
            {
              print(paste0("j=", j))
              
              # then, check all surrounding points for 3rd point k of contiguous defect
              for( n in m_i:m_f )
              {
                k <- helper[j,n]
                
                if( !( is.na( k ) ) && !( is.na( pdp[k] ) ) )
                {
                  if( pdp[k] <= 5 )
                  {
                    # check that 1st and 3rd points are not the same
                    if( i != k )
                    {
                      print(paste0("k=", k))
                      
                      return( TRUE )
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  return( FALSE )
}

#' Low-pressure Glaucoma Treatment Study Criteria (LOGTS)
#'
#' Analyse a VF for the presence of a VF defect based on the following criteria: at least 3 contiguous points depressed more than 8 decibels or 2 contiguous points depressed more than 10 decibels on total deviation plot
#' @param td total deviation probability map to be analyzed
#' @param helper helper map tht specifies the surrounding points of each point
#' @return BOOL Whether VF analysis resulted in a detection of a visual field defect
logts <- function( td, helper )
{
  print( "logts")
  
  # get index of first and last surround map points
  m_i <- which( colnames( helper ) == "topleft" )
  m_f <- which( colnames( helper ) == "botright" )
  
  # check cluster of 2 points
  # start by checking each point i in the td
  for( i in 1:length( td ) )
  {
    if( !( is.na( td[i] ) ) )
    {
      if( td[i] <= -10 )
      {
        print(paste0("i=", i))
        
        # then, check all surrounding points for 2nd point j of contiguous defect
        for( m in m_i:m_f )
        {
          j <- helper[i,m]
          
          if( !( is.na( j ) ) && !( is.na( td[j] ) ) )
          {
            if( td[j] <= -10 )
            {
              print(paste0("j=", j))
              
              return( TRUE )
            }
          }
        }
      }
    }
  }
  
  # check cluster of 3 points
  # start by checking each point i in the d
  for( i in 1:length( td ) )
  {
    if( !( is.na( td[i] ) ) )
    {
      if( td[i] <= -8 )
      {
        print(paste0("i=", i))
        
        # then, check all surrounding points for 2nd point j of contiguous defect
        for( m in m_i:m_f )
        {
          j <- helper[i,m]
          
          if( !( is.na( j ) ) && !( is.na( td[j] ) ) )
          {
            if( td[j] <= -8 )
            {
              print(paste0("j=", j))
              
              # then, check all surrounding points for 3rd point k of contiguous defect
              for( n in m_i:m_f )
              {
                k <- helper[j,n]
                
                if( !( is.na( k ) ) && !( is.na( td[k] ) ) )
                {
                  if( td[k] <= -8 )
                  {
                    # check that 1st and 3rd points are not the same
                    if( i != k )
                    {
                      print(paste0("k=", k))
                      
                      return( TRUE )
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  return( FALSE )
}

