#' @rdname vfevent
#' @author Herman Stubeda, Dalhousie University Medical School
#' @title Event criteria for detecting glaucomatous damage
#' @description Analyze a visual field for the presence of a defect based on
#'  the following criteria: Glaucoma Hemifield Test (GHT), Hodapp-Anderson-Parrish 2 (HAP2), 
#'  Foster, United Kingdom, Glaucoma Treatment Study (UKGTS), 
#'  and/or Low-pressure Glaucoma Treatment Study (LoGTS)
#' @details  
#'  The \code{vfdefect} function is a wrapper that contains all 5 event criteria. The user 
#'  can select which \code{criteria} to apply to the \code{vf} by passing \code{all} (default), 
#'  \code{ght}, \code{hap2}, \code{foster}, \code{ukgts}, or \code{logts}.
#'  
#'  The user can provide their own limits of normality \code{ghtlims} for GHT, or let the 
#'  function compute the limits by providing a vf object \code{vfctr} that contains normal visual fields.
#'  
#' @details Event criteria description 
#' \itemize{
#'   \item GHT: Decision tree focused on the comparison of 5 regions in the superior hemifield 
#'   with their mirror images in the inferior hemifield. Criteria commonly considered positive 
#'   if the output is "outside normal limits" [1]
#'   \item HAP2: GHT "Outside normal limits" \emph{OR} a cluster of 3 points depressed at P<0.05 level, 
#'   one of which at P<0.01 on the pattern deviation plot \emph{OR} PSD at P<5% [2]
#'   \item Foster: GHT "outside normal limits" \emph{AND} a cluster of 3 contiguous points depressed at 
#'   P<0.05 level on the pattern deviation plot [3]
#'   \item UKGTS: a cluster of 2 points depressed at P<0.01 level \emph{OR} a cluster of 3 points 
#'   depressed at P<0.05 level \emph{OR} a 10 dB difference across the nasal horizontal midline 
#'   at 2 or more adjacent points on the total deviation plot [4]
#'   \item LoGTS: a cluster of 3 points depressed at < -8 dB \emph{OR} a cluster of 2 
#'   points depressed at < -10 dB on the total deviation plot [5]
#' }
#' 
#' @param vf Visual fields to be analyzed as a standard vfobject
#' @param criteria select to use \code{all} (default) or one criteria:
#' \code{ght}, \code{hap2}, \code{foster}, \code{ukgts}, or \code{logts}
#' @param ghtlims a list containing limits of normality for the Glaucoma Hemifield 
#' Test (GHT); this list is returned by \code{vfdefect} function
#' @param vfctr control data, which is a vf object containing normative vfs
#' @param tpattern pattern of visual field test, either 24-2 or 30-2
#' @export
#' @return A list of 2 objects: 
#' \itemize{
#'  \item \code{$ght.limits} limits of normality required by the Glaucoma Hemifield Test (GHT) 
#'  \item \code{$results} results of detected glaucoma damage in \code{vf} parameter using the 
#'  user-selected \code{criteria} 
#' }
#' @examples
#'  # check visual fields for defect using all 5 criteria
#'  result <- vfdefect( vfpwgSunyiu24d2 )
#'  
#'  # check visual field for defect using the GHT criteria
#'  result <- vfdefect( vfpwgSunyiu24d2, criteria = "ght")
#'  
#'  # speed up computation when re-using the same ght limits of normality
#'  result1 <- vfdefect( vfpwgSunyiu24d2 )
#'  result2 <- vfdefect( vfpwgRetest24d2, ghtlims = result1$ght.limits )
#'  
#'  # use a custom normative dataset with the criteria to detect glaucoma defect:
#'  # 1.prepare the statistical environment with your normative dataset
#'  # 2.check visual fields for defect using all 5 criteria
#'  setnv( nvgenerate( vfctrSunyiu24d2 ) )
#'  result <- vfdefect( vfpwgSunyiu24d2, vfctr = vfctrSunyiu24d2 )
#' @references
#'  [1] Asman P, Heijl A. \emph{Glaucoma Hemifield Test: Automated Visual Field Evaluation.} 
#'  Arch Ophthalmol. 1992;110(6):812-819.
#'  
#'  [2] Hodapp E, Parrish RK, Anderson D. \emph{Clinical Decisions in Glaucoma.} 
#'  1st ed. Mosby, Inc.; 1993.
#'  
#'  [3] Foster PJ, Buhrmann R, Quigley HA, Johnson GJ. \emph{The definition and 
#'  classification of glaucoma in prevalence surveys.} Br J Ophthalmol. 2002;86(2):238-42.
#'  
#'  [4] Garway-Heath DF, Lascaratos G, Bunce C, Crabb DP, Russell RA, Shah A. 
#'  \emph{The United Kingdom glaucoma treatment study: A multicenter, randomized, 
#'  placebo-controlled clinical trial: Design and methodology.} Ophthalmology. 
#'  2013;120(1):68-76.
#'  
#'  [5] Krupin T, Liebmann JM, Greenfield DS, Rosenberg LF, Ritch R, Yang JW. \emph{The 
#'  Low-pressure Glaucoma Treatment Study (LoGTS): Study design and baseline characteristics 
#'  of enrolled patients.} Ophthalmology. 2005;112(3):376-385.

vfdefect <- function(vf , criteria = "all", ghtlims = NA, vfctr = visualFields::vfctrSunyiu24d2,  tpattern = "p24d2") {
  if( nrow( vf ) < 1 )
    stop( "vf is empty" )
  if( nrow( vfctr ) < 1 )
    stop( "control vf is empty" )
  if( is.list( ghtlims ) )
  {
    if( length( ghtlims$gh ) != 2 || length( ghtlims$sector[1,] ) != 6 || length( ghtlims$sector[,1] ) != 5 ) 
      stop( "ght limits object has improper structure" )
  }
  
  # compute td, pd maps
  td <- gettd( vf )
  pd <- getpd( td )

  if( tpattern == "p24d2" )
  {
    # defines the indices of surrounding points for each point on the 24-2 vf map
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
    
    # compute tdp, pdp, and extract field columns only
    range <- which( colnames( vf ) == "l1") : which( colnames( vf ) == "l54" )
    tdp <- gettdp( td )[range]
    pdp <- getpdp( pd )[range]
    td  <-         td  [range]
    pd  <-         pd  [range]
  }
  else if( tpattern == "p30d2" )
  {
    # defines the indices of surrounding points for each point on the 32-2 vf map
    SURR_MAP <- data.frame(
      "x" =          c(-9,-3,3,9,-15,-9,-3,3,9,15,-21,-15,-9,-3,3,9,15,21,-27,-21,-15,-9,-3,3,9,15,21,27,-27,-21,-15,-9,-3,3,9,15,21,27,-27,-21,-15,-9,-3,3,9,15,21,27,-27,-21,-15,-9,-3,3,9,15,21,27,-21,-15,-9,-3,3,9,15,21,-15,-9,-3,3,9,15,-9,-3,3,9),
      "y" =	         c(27,27,27,27,21,21,21,21,21,21,15,15,15,15,15,15,15,15,9,9,9,9,9,9,9,9,9,9,3,3,3,3,3,3,3,3,3,3,-3,-3,-3,-3,-3,-3,-3,-3,-3,-3,-9,-9,-9,-9,-9,-9,-9,-9,-9,-9,-15,-15,-15,-15,-15,-15,-15,-15,-21,-21,-21,-21,-21,-21,-27,-27,-27,-27),
      "topleft" =		 c(NA,NA,NA,NA,NA,NA,1,2,3,4,NA,NA,5,6,7,8,9,10,NA,NA,11,12,13,14,15,16,17,18,NA,19,20,21,22,23,24,25,26,27,NA,29,30,31,32,33,34,35,36,37,NA,39,40,41,42,43,44,45,46,47,49,50,51,52,53,54,55,56,59,60,61,62,63,64,67,68,69,70),
      "top"	=				 c(NA,NA,NA,NA,NA,1,2,3,4,NA,NA,5,6,7,8,9,10,NA,NA,11,12,13,14,15,16,17,18,NA,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,50,51,52,53,54,55,56,57,60,61,62,63,64,65,68,69,70,71),
      "topright"	=	 c(NA,NA,NA,NA,1,2,3,4,NA,NA,5,6,7,8,9,10,NA,NA,11,12,13,14,15,16,17,18,NA,NA,20,21,22,23,24,25,26,27,28,NA,30,31,32,33,34,35,36,37,38,NA,40,41,42,43,44,45,46,47,48,NA,51,52,53,54,55,56,57,58,61,62,63,64,65,66,69,70,71,72),
      "left"	=	     c(NA,1,2,3,NA,5,6,7,8,9,NA,11,12,13,14,15,16,17,NA,19,20,21,22,23,24,25,26,27,NA,29,30,31,32,33,34,35,36,37,NA,39,40,41,42,43,44,45,46,47,NA,49,50,51,52,53,54,55,56,57,NA,59,60,61,62,63,64,65,NA,67,68,69,70,71,NA,73,74,75),
      "right"	=      c(2,3,4,NA,6,7,8,9,10,NA,12,13,14,15,16,17,18,NA,20,21,22,23,24,25,26,27,28,NA,30,31,32,33,34,35,36,37,38,NA,40,41,42,43,44,45,46,47,48,NA,50,51,52,53,54,55,56,57,58,NA,60,61,62,63,64,65,66,NA,68,69,70,71,72,NA,74,75,76,NA),
      "botleft"	=    c(5,6,7,8,11,12,13,14,15,16,19,20,21,22,23,24,25,26,NA,29,30,31,32,33,34,35,36,37,NA,39,40,41,42,43,44,45,46,47,NA,49,50,51,52,53,54,55,56,57,NA,NA,59,60,61,62,63,64,65,66,NA,NA,67,68,69,70,71,72,NA,NA,73,74,75,76,NA,NA,NA,NA),
      "bot"	=        c(6,7,8,9,12,13,14,15,16,17,20,21,22,23,24,25,26,27,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,NA,59,60,61,62,63,64,65,66,NA,NA,67,68,69,70,71,72,NA,NA,73,74,75,76,NA,NA,NA,NA,NA),			
      "botright" =   c(7,8,9,10,13,14,15,16,17,18,21,22,23,24,25,26,27,28,30,31,32,33,34,35,36,37,38,NA,40,41,42,43,44,45,46,47,48,NA,50,51,52,53,54,55,56,57,58,NA,59,60,61,62,63,64,65,66,NA,NA,67,68,69,70,71,72,NA,NA,73,74,75,76,NA,NA,NA,NA,NA,NA),
      "refl" =       c(73,74,75,76,67,68,69,70,71,72,59,60,61,62,63,64,65,66,49,50,51,52,53,54,55,56,57,58,39,40,41,42,43,44,45,46,47,48,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),
      "ght.sector" = c(NA,NA,NA,NA,NA,1,1,2,2,NA,NA,1,1,1,1,2,2,NA,NA,3,3,4,4,4,4,NA,NA,NA,3,3,3,5,5,5,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA)
    )		
    
    # compute tdp, pdp, and extract relevant columns
    range <- which( colnames( vf ) == "l1" ) : which( colnames( vf ) == "l76" )
    tdp <- gettdp( td )[range]
    pdp <- getpdp( pd )[range]
    td  <-         td  [range]
    pd  <-         pd  [range]
  }
  else
  {
    stop( "invalid tperimtery test type: pass only p24d2 or p32d2" )
  }

  if( !is.list( ghtlims ) )
  {
    #--------------------------------------------------------------------------------------------------------------------------------------   
    # Calculate limits of normality as described in Asman & Heijl, Arch Opthalmol, 1992
    # 1. pass about 100+ vfs of normal eyes, preferably with 1 unique eye per patient, within the vf_ctrl parameter
    # 2. for each vf in vf_ctrl, calculate gh, sector sums, and sector up-down differences
    # 3. use the bootstrap method to calculate limits of normality for ght, sector sums, and sector up-down differences
    #    a. select 500 random samples from vf_ctrl, making sure each sample has 1 unique eye per patient
    #    b. for each sample, calculate the confidence intervals: 99% for gh, 99% for sector sums, 99% and 97% for sector up-down differences
    #    c. compute means of limits from the 500 samples
  
    # define constants
    CTRL_VF_NUM     <- nrow(vfctr)
    SUMS_NUM        <- 3
    GHT_SECTORS_NUM <- 5
    SECTORS_LIM_NUM <- 6
    GH_LIM_NUM      <- 2
    SAMPLE_NUM      <- 500
    
    # compute td, pd, pdp maps, and extract relevant rows
    tdctr  <- gettd ( vfctr )
    pdctr  <- getpd ( tdctr )
    pdpctr <- getpdp( pdctr )[range]
    pdctr  <-         pdctr  [range]
    
    # define data structures for values of gh, sector sums, and sector up-down differences
    sector.sums <- array(NA, c( GHT_SECTORS_NUM, SUMS_NUM, CTRL_VF_NUM ) )
    gh          <- array(NA, c(               1,        1, CTRL_VF_NUM ) )
    
    colnames( sector.sums ) <- c( "sum.sup", "sum.inf", "up.down" )
    colnames(          gh ) <- c( "gh" )
    
    # compute gh, sector sums, and sector up-down differences
    for( v in 1:CTRL_VF_NUM )
    {
      scores <- pdpctr[v,]
      scores[1,] <- 0
      
      for( i in 1:length( pdpctr[v,] ) )
      {
        if( !( is.na( pdpctr[v,i] ) ) )
        {
          if( pdpctr[v,i] > 5 )
            scores[i] <- 0
          else if( pdpctr[v,i] > 2 )
            scores[i] <- 2
          else if( pdpctr[v,i] > 1 )
            scores[i] <- 5
          else
            scores[i] <- 10 * abs( pdctr[v,i] / getnv()$luts$pd["1%",i] )  
        }
      }
      for( sector in 1:GHT_SECTORS_NUM )
      {
        sector.sums[ sector,"sum.sup",v ] <- sum( scores[ which( SURR_MAP$ght.sector == sector ) ] )
        sector.sums[ sector,"sum.inf",v ] <- sum( scores[ SURR_MAP[ which( SURR_MAP$ght.sector == sector ), "refl" ] ] )
        sector.sums[ sector,"up.down",v ] <- sector.sums[ sector,"sum.sup",v ] - sector.sums[ sector,"sum.inf",v ]
      }
      
      gh[ 1,"gh", ] <- getgh( tdctr )
    }
    
    # define data structure for limits of gh, sector sums, and sector up-down differences
    sector.lims <- data.frame( array( NA, c( GHT_SECTORS_NUM, SECTORS_LIM_NUM ) ) )
    gh.lims     <- data.frame( array( NA, c(               1,      GH_LIM_NUM ) ) )
    
    colnames( sector.lims ) <- c("sum.sup99.5", "sum.inf99.5", "up.down0.5", "up.down99.5", "up.down1.5", "up.down98.5" )
    colnames( gh.lims )     <- c( "gh0.5", "gh99.5" )
    
    # compute limits of gh, sector sums, and sector up-down differences
    for( sector in 1:GHT_SECTORS_NUM )
    {
      sector.lims[ sector,"sum.sup99.5"] <- boot( data=sector.sums[ sector,"sum.sup", ], statistic=lim, R=SAMPLE_NUM, l=0.995 )$t0
      sector.lims[ sector,"sum.inf99.5"] <- boot( data=sector.sums[ sector,"sum.inf", ], statistic=lim, R=SAMPLE_NUM, l=0.995 )$t0
      sector.lims[ sector,"up.down99.5"] <- boot( data=sector.sums[ sector,"up.down", ], statistic=lim, R=SAMPLE_NUM, l=0.995 )$t0
      sector.lims[ sector,"up.down0.5" ] <- boot( data=sector.sums[ sector,"up.down", ], statistic=lim, R=SAMPLE_NUM, l=0.005 )$t0
      sector.lims[ sector,"up.down98.5"] <- boot( data=sector.sums[ sector,"up.down", ], statistic=lim, R=SAMPLE_NUM, l=0.985 )$t0
      sector.lims[ sector,"up.down1.5" ] <- boot( data=sector.sums[ sector,"up.down", ], statistic=lim, R=SAMPLE_NUM, l=0.015 )$t0
    }
    gh.lims$gh99.5 <- boot( data=gh[ 1,"gh", ], statistic=lim, R=SAMPLE_NUM, l=0.995)$t0
    gh.lims$gh0.5  <- boot( data=gh[ 1,"gh", ], statistic=lim, R=SAMPLE_NUM, l=0.005)$t0
    
    # assign data structure with all computed limits
    ght.lims <- list( "gh" = gh.lims, "sector" = sector.lims )
  }
  else
  {
    ght.lims <- ghtlims
  }
  
  #---------------------------------------------------------------------------------------------------------------------------------------
  #Apply vf criteria to detect glaucoma defect

  # implement the results data frame
  if( criteria == "all")
  {
    result <- data.frame( array( FALSE, c( nrow( vf ), 5 ) ) )
    colnames( result ) <- c( "ght", "hap2", "foster", "ukgts", "logts" )
  }
  else
  {
    result <- data.frame( array( FALSE, c( nrow( vf ), 1 ) ) )
    colnames( result ) <- criteria
  }
  rownames( result ) <- rownames( vf )
  
  # analyze vf for defect
  for( i in 1:nrow(vf) )
  {
    #print(i)
    if( criteria == "all" )
    {
      ght <- ght( vf[i,], td[i,], pd[i,], pdp[i,], ght.lims, SURR_MAP )
      
      result[i,] <- c( ght, 
                       hap2  ( vf [i,], pd [i,], pdp[i,], ght, SURR_MAP ), 
                       foster(          pd [i,], pdp[i,], ght, SURR_MAP ), 
                       ukgts ( vf [i,], td [i,], tdp[i,], SURR_MAP , tpattern ), 
                       logts (          td [i,],          SURR_MAP ) 
      )
    }
    
    else if( criteria == "ght" )
      result[i,] <- ght( vf[i,], td[i,], pd[i,], pdp[i,], ght.lims, SURR_MAP )    
    
    else if( criteria == "hap2" )
    {
      ght <- ght( vf[i,], td[i,], pd[i,], pdp[i,], ght.lims, SURR_MAP )
      
      result[i,] <- hap2( vf[i,], pdp[i,], ght, SURR_MAP)
    }
    
    else if( criteria == "foster" )
    {
      ght <- ght( vf[i,], td[i,], pd[i,], pdp[i,], ght.lims, SURR_MAP )
      
      result[i,] <- foster( pdp[i,], ght, SURR_MAP )
    } 
    
    else if( criteria == "ukgts" )
      result[i,] <- ukgts( vf[i,], td[i,], tdp[i,], SURR_MAP, tpattern )
    
    else if( criteria == "logts" )
      result[i,] <- logts( td[i,], SURR_MAP )
    
    else
      stop( "invalid criteria: select all (default), hap2, ukgts, ght, foster, or logts" )
  }
  
  return( list( "ght.limits" = ght.lims, "results" = result) )
}

####################
# INTERNAL FUNCTIONS
####################
# description: Compute the limit of probability
# param: data a one-dimensional vector of data from which to compute the quantile limit
# param: indices a variable that is required by the boot function, which will be used to sample the data parameter
# param: l quantile
#' @noRd
lim <- function(data, indices, l)
  return( quantile( data[indices], p=l, names=FALSE ) )

# description: visual field defect criteria for detecting glaucomatous damage
# param vf Visual fields to be analyzed as a standard vfobject
# param td total deviation map to be analyzed
# param tdp total deviation probability map to be analyzed
# param pd pattern deviation map to be analyzed
# param pdp pattern deviation probability map to be analyzed
# param lims limits of normality required for the GHT algorithm
# param helper helper map tht specifies the surrounding points of each point
# tpattern pattern of visual field test, either 24-2 or 30-2
#' @noRd
ght <- function(vf, td, pd, pdp, lims, helper)
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
      score[i] <- 10 * abs( pd[i] / getnv()$luts$pd["1%",i] )  
    
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

#' @noRd
hap2 <- function(vf, pd, pdp, ght, helper)
{
  #print( "hap2")
  
  # check GHT result of vf
  if( "Outside normal limits" == ght )
    return( TRUE )
  
  # check PSD result of vf
  psd <- getglp( getgl( vf ) )["psd"]
  if( !is.na(psd) )
  {
    if( psd <= 5 )
      return( TRUE )
  }
  
  # if no pd map exists for a given vf, then depression in td must be severe enough to indicate a hypothetical cluster on pdp
  if( length( which( is.na( pd ) ) ) == 54 )
  {
    return( TRUE )
  }
  
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
        #print(paste0("i=", i))
        
        # then, check all surrounding points for 2nd point j of contiguous defect
        for( m in m_i:m_f )
        {
          j <- helper[i,m]
          
          if( !( is.na( j ) ) && !( is.na( pdp[j] ) ) )
          {
            if( pdp[j] <= 5 )
            {
              #print(paste0("j=", j))
              
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
                        #print(paste0("k=", k))
                        
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

#' @noRd
foster <- function(pd, pdp, ght, helper)
{
  #print( "foster" )
  
  # obtain GHT result of vf
  if( "Outside normal limits" != ght )
    return( FALSE )
  
  # if no pd map exists for a given vf, then depression in td must be severe enough to indicate a hypothetical cluster on pdp
  if( length( which( is.na( pd ) ) ) == 54 )
  {
    return( TRUE )
  }
  
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
        #print(paste0("i=", i))
        
        # then, check all surrounding points for 2nd point j of contiguous defect
        for( m in m_i:m_f )
        {
          j <- helper[i,m]
          
          if( !( is.na( j ) ) && !( is.na( pdp[j] ) ) )
          {
            if( pdp[j] <= 5 )
            {
              #print(paste0("j=", j))
              
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
                      #print(paste0("k=", k))
                      
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

#' @noRd
ukgts <- function(vf, td, tdp, helper, tpattern)
{
  #print( "ukgts")
  
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
        #print(paste0("i=", i))
        
        # then, check all surrounding points for 2nd point j of contiguous defect
        for( m in m_i:m_f )
        {
          j <- helper[i,m]
          
          if( !( is.na( j ) ) && !( is.na( tdp[j] ) ) )
          {
            if( tdp[j] <= 1 )
            {
              #print(paste0("j=", j))
              
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
        #print(paste0("i=", i))
        
        # then, check all surrounding points for 2nd point j of contiguous defect
        for( m in m_i:m_f )
        {
          j <- helper[i,m]
          
          if( !( is.na( j ) ) && !( is.na( tdp[j] ) ) )
          {
            if( tdp[j] <= 5 )
            {
              #print(paste0("j=", j))
              
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
                      #print(paste0("k=", k))
                      
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
  if( tpattern == "p24d2" )
  {
    ir_i <- which( colnames( td ) == "l19" )
    ir_f <- which( colnames( td ) == "l23" )
  }
  else
  {
    ir_i <- which( colnames( td ) == "l29" )
    ir_f <- which( colnames( td ) == "l33" )    
  }
  
  
  # check reflection points
  # start by checking each point i in the td
  for( i in ir_i:ir_f )
  {
    ir <- helper[i,"refl"]
    
    if( !( is.na( td[i] ) ) && !( is.na( td[ir] ) ) )
    {
      if( abs( td[i] - td[ir] ) >= 10 )
      {
        #print(paste0("i=", i, " ir=", ir))
        
        j <- helper[i,"right"]
        
        jr <- helper[j,"refl"]
        
        if( !( is.na( j ) ) && !( is.na( jr ) ) )
        {
          if( abs( td[j] - td[jr] ) >= 10 )
          {
            if( ( i != j ) && ( ir != j ) )
            {
              #print(paste0("j=", j, " jr=", jr))
              
              return( TRUE )
            }
          }
        }
      }
    }
  }
  return(FALSE)
}

#' @noRd
logts <- function(td, helper)
{
  #print( "logts")
  
  # get index of first and last surround map points
  m_i <- which( colnames( helper ) == "topleft" )
  m_f <- which( colnames( helper ) == "botright" )
  
  # check cluster of 2 points
  # start by checking each point i in the td
  for( i in 1:length( td ) )
  {
    if( !( is.na( td[i] ) ) )
    {
      if( td[i] < -10 )
      {
        #print(paste0("i=", i))
        
        # then, check all surrounding points for 2nd point j of contiguous defect
        for( m in m_i:m_f )
        {
          j <- helper[i,m]
          
          if( !( is.na( j ) ) && !( is.na( td[j] ) ) )
          {
            if( td[j] < -10 )
            {
              #print(paste0("j=", j))
              
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
      if( td[i] < -8 )
      {
        #print(paste0("i=", i))
        
        # then, check all surrounding points for 2nd point j of contiguous defect
        for( m in m_i:m_f )
        {
          j <- helper[i,m]
          
          if( !( is.na( j ) ) && !( is.na( td[j] ) ) )
          {
            if( td[j] < -8 )
            {
              #print(paste0("j=", j))
              
              # then, check all surrounding points for 3rd point k of contiguous defect
              for( n in m_i:m_f )
              {
                k <- helper[j,n]
                
                if( !( is.na( k ) ) && !( is.na( td[k] ) ) )
                {
                  if( td[k] < -8 )
                  {
                    # check that 1st and 3rd points are not the same
                    if( i != k )
                    {
                      #print(paste0("k=", k))
                      
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