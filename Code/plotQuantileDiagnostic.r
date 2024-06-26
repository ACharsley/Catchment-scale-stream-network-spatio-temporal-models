###############################################################################################################################
##
##  plotQuantileDiagnostic function
##
###############################################################################################################################

######## Define the "plotQuantileDiagnostic" function
plotQuantileDiagnostic <- function( TmbData, Report, DateFile = paste0( getwd(), "/" ), 
	save_dir = paste0( DateFile, "/QQ_Fn/" ), FileName_PP = "Posterior_Predictive",
	FileName_Phist = "Posterior_Predictive-Histogram", FileName_QQ = "Q-Q_plot", FileName_Qhist = "Q-Q_hist" ) {

	#### Retrieve data based on model type
    	if( "n_e" %in% names( TmbData ) ) {
		## VAST Versions 3.0.0, 4.0.0: names n_e, e_i and ObsModel_ez are added to support for error distributions
        	n_e <- TmbData$n_e
        	e_i <- TmbData$e_i
        	ObsModel_ez <- TmbData$ObsModel_ez
        	sigmaM <- Report$SigmaM
    	} else if( "n_c" %in% names( TmbData ) ) {
        	## VAST Version < 3.0.0: names n_c, c_i and ObsModel are used to group by categories
        	n_e <- TmbData$n_c
        	e_i <- TmbData$c_i
        	ObsModel_ez <- rep( 1, n_e ) %o% TmbData$ObsModel
        	sigmaM <- Report$SigmaM
        	if( is.vector( sigmaM ) ) 
            	## VAST Versions 1.0.0 and 1.1.0
            	sigmaM <- rep( 1, n_e ) %o% Report$SigmaM
    	} else {
        	## SpatialDeltaGLMM
        	n_e <- 1
        	e_i <- rep( 0, TmbData$n_i )
        	ObsModel_ez <- matrix( TmbData$ObsModel, nrow = 1 )
        	sigmaM <- matrix( Report$SigmaM, nrow = 1 )
    	}

	#### Check data
    	if( nlevels( as.factor( e_i )) != n_e ) stop( "Error in e_i: nlevels does not agree with n_e" )
    	if( nrow( as.matrix( ObsModel_ez ) ) != n_e ) stop( "Error in ObsModel_ez: nrow does not agree with n_e" )
    	if( nrow( as.matrix( sigmaM ) ) != n_e ) stop( "Error in sigmaM: nrow does not agree with n_e" )
    
	#### Check save directory
    	dir.create( save_dir, recursive = TRUE, showWarnings = FALSE )
    	if( !dir.exists( save_dir ) ) stop( paste0( "Wrong directory, cannot save plots: ", save_dir ) )
    
	#### Return list
    	Return <- vector( "list", length = n_e )
  
	#### Define an utility function
    	pow = function( a, b ) a^b

    	#### Loop through each group (plot functions remain unchanged from previous version)
    	for( i_e in 1 : n_e ){

      	## Generate plot names
      	if( !is.null( save_dir ) ) {
        		if( !is.null( FileName_PP ) ) save_PP = paste0( save_dir, "/", FileName_PP, "-", i_e, ".jpg" )
        		if( !is.null( FileName_Phist ) ) save_Phist = paste0( save_dir, "/", FileName_Phist, "-", i_e, ".jpg" )
        		if( !is.null( FileName_QQ ) ) save_QQ = paste0( save_dir, "/", FileName_QQ, "-", i_e, ".jpg" )
        		if( !is.null( FileName_Qhist ) ) save_Qhist = paste0( save_dir, "/", FileName_Qhist, "-", i_e, ".jpg" )
      	}

      	## Find where b_i > 0 within category i_e
     		Which = which( TmbData$b_i > 0 & e_i == ( i_e - 1 ) )
      	Q = rep( NA, length( Which ) ) # Vector to track quantiles for each observation
      	y = array( NA, dim = c( length( Which ), 1000 ) ) # Matrix to store samples
      	pred_y = var_y = rep( NA, length( Which ) ) # Vector to track quantiles for each observation

      	## Calculate pred_y
      	## We cannot use R2_i anymore because interpretation changed around March 9, 2017 (due to area-swept change 
		## in Poisson-process and Tweedie functions). However, we can use P2_i, which has a stable definition over 
		## time (as a linear predictor)
      	if( length( ObsModel_ez[i_e,] ) >= 2 && ObsModel_ez[i_e,2] == 2 ){
        	Return[[i_e]] = list( "type" = ObsModel_ez[i_e,], message = "QQ not set up for Tweedie distribution" )
       	next
	}

      if( !(ObsModel_ez[i_e,1] %in% c(1,2)) ){
        Return[[i_e]] = list("type"=ObsModel_ez[i_e,], message="QQ not working except for when using a Gamma or Lognormal distribution")
        next
      }
      if( length( ObsModel_ez[i_e,] ) == 1 || ObsModel_ez[i_e,2] %in% c( 0, 3 ) ) {
        	for( ObsI in 1 : length( Which ) ){
          		pred_y[ObsI] = TmbData$a_i[Which[ObsI]] * exp( Report$P2_i[Which[ObsI]] )
        	}
      }
      if( length( ObsModel_ez[i_e,] )>= 2 && ObsModel_ez[i_e,2] %in% c( 1, 4 ) ) {
		for( ObsI in 1 : length( Which ) ) {
          		if( sigmaM[e_i[Which[ObsI]]+1,3] != 1 ) {
				stop( "`QQ_Fn` will not work with the Poisson-link delta model across all VAST versions given values for turned-off parameters" )
			}
          		R1_i = 1 - exp( -1 * sigmaM[e_i[Which[ObsI]]+1,3] * TmbData$a_i[Which[ObsI]] * 
				exp( Report$P1_i[Which[ObsI]] ) )
          		pred_y[ObsI] = TmbData$a_i[Which[ObsI]] * exp( Report$P1_i[Which[ObsI]] ) / 
				R1_i * exp( Report$P2_i[Which[ObsI]] ) ;
        	}
	}

	#### Simulate quantiles for different distributions: Loop through observations
      for( ObsI in 1 : length( Which ) ) {
		if( ObsModel_ez[i_e,1] == 1 ) {
			y[ObsI,] = rlnorm( n = ncol( y ), meanlog = log( pred_y[ObsI] ) - pow( sigmaM[i_e,1], 2 ) / 2, 
				sdlog = sigmaM[i_e,1] ) ## Plotting in log-space
          		Q[ObsI] = plnorm( q = TmbData$b_i[Which[ObsI]], meanlog = log( pred_y[ObsI] ) - pow( sigmaM[i_e,1], 2 ) / 2, 
				sdlog = sigmaM[i_e,1] )
        	}
        	if( ObsModel_ez[i_e,1] == 2 ) {
          		b = pow( sigmaM[i_e, 1], 2 ) * pred_y[ObsI] ;
          		y[ObsI,] = rgamma( n = ncol( y ), shape = 1 / pow( sigmaM[i_e,1] , 2 ), scale = b )
          		Q[ObsI] = pgamma( q = TmbData$b_i[Which[ObsI]], shape = 1 / pow( sigmaM[i_e,1], 2 ), scale = b )
        	}
	}

      #### Make plot while calculating posterior predictives
      if( !is.null( FileName_PP ) & !is.null( save_dir ) ) jpeg( save_PP, width = 10, height = 3, res = 200, units = "in" )
      par( mar = c( 2, 2, 2, 0 ), mgp = c( 1.25, 0.25, 0 ), tck = -0.02 )
      plot( TmbData$b_i[Which], ylab = "", xlab = "", log = "y", main = "", col = "blue")

      #### Add results to plot: Loop through observations
      for( ObsI in 1 : length( Which ) ) {
		var_y[ObsI] = var( y[ObsI,] )
        	Quantiles = quantile( y[ObsI,], prob = c( 0.025 ,0.25, 0.75, 0.975 ) )
        	lines( x = c( ObsI, ObsI ), y = Quantiles[2:3], lwd = 2 )
        	lines( x = c( ObsI, ObsI ), y = Quantiles[c( 1, 4 )], lwd = 1, lty = "dotted" )
        	if( TmbData$b_i[Which[ObsI]] > max(Quantiles) | TmbData$b_i[Which[ObsI]] < min( Quantiles ) ) {
          		points( x = ObsI, y = TmbData$b_i[Which[ObsI]], pch = 4, col = "red", cex = 2 )
        	}
      }
      if( !is.null( FileName_PP ) & !is.null( save_dir ) ) dev.off()

      #### Produce a Q-Q plot
      if( !is.null( FileName_QQ ) & !is.null( save_dir ) ) jpeg( save_QQ, width = 4, height = 4, res = 200, units = "in" )
      par( mfrow = c( 1, 1 ), mar = c( 2, 2, 2, 0 ), mgp = c( 1.25, 0.25, 0 ), tck = -0.02 )
      Qtemp = na.omit( Q )
      Order = order( Qtemp )
      plot( x = seq( 0, 1, length = length( Order ) ), y = Qtemp[Order], main = "Q-Q plot", xlab = "Uniform", 
		ylab = "Empirical", type = "l", lwd = 3 )
      abline( a = 0, b = 1 )
      if( !is.null( FileName_QQ ) & !is.null( save_dir ) ) dev.off()

      #### Plot aggregate predictive distribution
      if( !is.null( FileName_Phist ) & !is.null( save_dir )) jpeg( save_Phist, width = 4, height = 4, res = 200, units = "in" )
      par( mfrow = c( 1, 1 ), mar = c( 2, 2, 2, 0 ), mgp = c( 1.25, 0.25, 0 ), tck = -0.02 )
      hist( log(y), main = "Aggregate predictive distribution", xlab = "log( Obs )", ylab = "Density" )
      if( !is.null( FileName_Phist ) & !is.null( save_dir ) ) dev.off()

      #### Produce a quantile histogram
      if( !is.null( FileName_Qhist ) & !is.null( save_dir ) ) jpeg( save_Qhist, width = 4, height = 4, res = 200, units = "in")
      par( mfrow = c( 1, 1 ), mar = c( 2, 2, 2, 0 ), mgp = c( 1.25, 0.25, 0 ), tck = -0.02 )
      hist( na.omit( Q ), main = "Quantile_histogram", xlab = "Quantile", ylab = "Number" )
      if( !is.null( FileName_Qhist ) & !is.null( save_dir )) dev.off()

      #### Return stuff
      Return[[i_e]] = list( "type" = ObsModel_ez[i_e,], "Q" = Q, "var_y" = var_y, "pred_y" = pred_y )
	}
    
	if( length( Return ) == 1 ) Return <- Return[[1]] ## Single species model
    	return( Return )
}
