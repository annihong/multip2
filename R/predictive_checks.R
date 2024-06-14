#' Calculate descriptive statistics for a network
#'
#' This function calculates the edge density, reciprocity, and transitivity of a network.
#'
#' @param network An igraph object representing the network.
#' @return A numeric vector of length 3. The first element is the edge density, the second element is the reciprocity, and the third element is the transitivity.
#' @importFrom igraph edge_density reciprocity transitivity
#' @examples
#' # TODO: Add examples of usage
#' @export
descriptive_stats_prior <- function(X){
    network <- igraph::graph_from_adjacency_matrix(X)
    mu <- igraph::edge_density(network)
    rho <- igraph::reciprocity(network)
    transitivity <- igraph::transitivity(network)
    res <- c(mu, rho, transitivity)
    return(res)
}

#' Perform checks on simulated networks
#'
#' This function takes a list of simulated networks, calculates descriptive statistics for each network,
#' and returns a ggplot2 object of histograms of these statistics.
#'
#' @param simulated_networks A list of simulated networks. Each element of the list should be an list of adjacency matrix (a multiplex network).
#' @param descriptive_stats A function that calculates descriptive statistics for a network. This function should take an igraph object as input and return a numeric vector of statistics.
#' @param stats_lab A character vector of labels for the statistics calculated by `descriptive_stats`. The default is `c("density", "reciprocity", "transitivity")`.
#' @return A ggplot2 object. This plot shows histograms of the calculated statistics for the simulated networks.
#' @import ggplot2
#' @importFrom igraph graph_from_adjacency_matrix
#' @importFrom reshape2 melt
#' @examples
#' # TODO: Add examples of usage
#' @export
prior_simulated_network_checks <- function(simulated_networks, descriptive_stats_func = descriptive_stats_prior, stats_lab = c("density", "reciprocity", "transitivity")) {
    plots <- list()
    for (layer_lab in names(simulated_networks)) {
        uniplex_nets <- simulated_networks[[layer_lab]]
        final_res <- list()
        for (uniplex_net in uniplex_nets) {
            res <- descriptive_stats_func(uniplex_net)
            final_res <- append(final_res, list(res))
        }
        final_res <- do.call(rbind, final_res)
        final_res[is.na(final_res)] <- 0
        colnames(final_res) <- stats_lab
        long <- as.data.frame(reshape2::melt(final_res)[,2:3])
        colnames(long) <- c("key", "val")
        dev.new()
        p <- ggplot(long, aes(x = val)) +
            geom_histogram() +
            facet_grid(.~key) +
            aes(y=after_stat(count)/sum(after_stat(count))) +
            #ylim(0, 0.75) +
            labs(title = layer_lab, x = "", y = "")+
            theme_bw() 
        plots <- append(plots, list(p))
        }
        do.call(gridExtra::grid.arrange, plots)
       
    
}

degree_distribution <- function(x, by_axis, levels, cumulative) {
    a <- apply(x, by_axis, sum)
	if (cumulative) {
		dd <- sapply(levels, function(i){ sum(a<=i) })
	}
	else {
		dd <- sapply(levels, function(i){ sum(a==i) })
	}
	names(dd) <- as.character(levels)
	return(dd)
}

#' Calculate the outdegree distribution of a single network
#'
#' This function calculates the outdegree distribution of a single network.
#' The outdegree distribution is a measure of the number of outgoing edges from each node in the network.
#'
#' @param single_net The single network for which to calculate the outdegree distribution.
#' @param levels The levels of outdegree to consider. Default is 0:8.
#' @param cumulative Logical indicating whether to calculate the cumulative distribution. Default is TRUE.
#'
#' @return The outdegree distribution of the single network.
#' @export
Outdegree_distribution <- function(single_net, levels=0:8, cumulative=TRUE) {
    return(degree_distribution(single_net, 1, levels, cumulative))
}
#' Calculate the indegree distribution of a single network
#'
#' This function calculates the indegree distribution of a single network.
#' The indegree distribution is a measure of the number of outgoing edges from each node in the network.
#'
#' @param single_net The single network for which to calculate the indegree distribution.
#' @param levels The levels of indegree to consider. Default is 0:8.
#' @param cumulative Logical indicating whether to calculate the cumulative distribution. Default is TRUE.
#'
#' @return The indegree distribution of the single network.
#' @export
Indegree_distribution <- function(single_net, levels=0:8, cumulative=TRUE) {
    return(degree_distribution(single_net, 2, levels, cumulative))
}

#' Triad Census
#'
#' This function calculates the triad census for a given network.
#'
#' @param single_net The adjacency matrix of the network.
#' @param levels The levels of triad types to include in the census.
#'   Default is 1:16, which includes all triad types.
#'
#' @return A vector containing the counts of each triad type specified by the levels parameter.
#' @export
Triad_census <- function (single_net, levels = 1:16) {
  # get matrix and prepare data
  mat <- single_net
  N <- nrow(mat)
  # matrix with reciprocal ties
  matReciprocal <- mat + t(mat)
  # matrix with direction information for triad
  matDirected <- mat - t(mat)
  matDirected[matDirected == -1] <- 2
  matDirected[matReciprocal == 2] <- 3
  matDirected <- matDirected + 1
  # reciproal matrix with ties from higher to lower IDs
  matHigher <- matReciprocal
  matHigher[lower.tri(matHigher)] <- 0
  # neighbors lookup
  neighbors<- apply(matReciprocal, 1, function(x) which(x > 0))
  # neighbors with lower ids
  neighborsHigher <- apply(matHigher, 1, function(x) which(x > 0))

  # lookup table for 64 triad types
  # i->j, j->k, i->k
  # 1: empty, 2: forward, 3: backward, 4: reciprocal
  # order as in vector tc
  lookup <- array(NA, dim = rep(4,3))
  lookup[1,1,1] <- 1
  lookup[2,1,1] <- lookup[1,2,1] <- lookup[1,1,2] <- lookup[3,1,1] <- lookup[1,3,1] <- lookup[1,1,3] <- 2
  lookup[4,1,1] <- lookup[1,4,1] <- lookup[1,1,4] <- 3
  lookup[2,1,2] <- lookup[3,2,1] <- lookup[1,3,3] <- 4
  lookup[2,3,1] <- lookup[3,1,3] <- lookup[1,2,2] <- 5
  lookup[2,2,1] <- lookup[3,3,1] <- lookup[2,1,3] <- lookup[3,1,2] <- lookup[1,2,3] <- lookup[1,3,2] <- 6
  lookup[4,3,1] <- lookup[4,1,3] <- lookup[2,4,1] <- lookup[1,4,2] <- lookup[3,1,4] <- lookup[1,2,4] <- 7
  lookup[4,2,1] <- lookup[4,1,2] <- lookup[3,4,1] <- lookup[1,4,3] <- lookup[2,1,4] <- lookup[1,3,4] <- 8
  lookup[2,2,2] <- lookup[2,3,3] <- lookup[2,3,2] <- lookup[3,3,3] <- lookup[3,2,2] <- lookup[3,2,3] <- 9
  lookup[2,2,3] <- lookup[3,3,2] <- 10  # 3-cycle
  lookup[4,4,1] <- lookup[4,1,4] <- lookup[1,4,4] <- 11
  lookup[2,4,2] <- lookup[3,2,4] <- lookup[4,3,3] <- 12
  lookup[2,3,4] <- lookup[3,4,3] <- lookup[4,2,2] <- 13
  lookup[2,2,4] <- lookup[3,3,4] <- lookup[2,4,3] <- lookup[3,4,2] <- lookup[4,2,3] <- lookup[4,3,2]<- 14
  lookup[2,4,4] <- lookup[4,2,4] <- lookup[4,4,2] <- lookup[3,4,4] <- lookup[4,3,4] <- lookup[4,4,3] <- 15
  lookup[4,4,4] <- 16

  # initialize triad census
  tc <- c("003"  = 0,
          "012"  = 0,
          "102"  = 0,
          "021D" = 0,
          "021U" = 0,
          "021C" = 0,
          "111D" = 0,
          "111U" = 0,
          "030T" = 0,
          "030C" = 0,
          "201"  = 0,
          "120D" = 0,
          "120U" = 0,
          "120C" = 0,
          "210"  = 0,
          "300"  = 0)

  # iterate through all non-empty dyads (from lower to higher ID)
  if (length(neighborsHigher) > 0){ # else mat is the zero matrix
	for(ii in 1:N){
		for(j in neighborsHigher[[ii]]){
      # set of nodes that are linked to ii and j
		third <- setdiff( union(neighbors[[ii]], neighbors[[j]]),
                    c(ii, j) )
		# store triads with just one tie
		triadType <- ifelse(matReciprocal[ii,j] == 2, 3, 2)
		tc[triadType] <- tc[triadType] + N - length(third) - 2
		for (k in third){
        # only store triads once
			if(j < k || ( ii < k && k < j && !(k %in% neighbors[[ii]]) ) ){
				t1 <- matDirected[ii,j]
				t2 <- matDirected[j,k]
				t3 <- matDirected[ii,k]
				triadType <- lookup[t1, t2, t3]
				tc[triadType] <- tc[triadType] + 1
				}
			}
		}
	}
  }
  # assign residual to empty triad count
  tc[1] <- 1/6 * N*(N-1)*(N-2) - sum(tc[2:16])
  return(tc[levels])
}

#' multiplex_gof_baseline
#'
#' This function compares the observed and simulated basic multiplex statistics using boxplots and scatter plots.
#'
#' @param dep_net The observed multiplex network.
#' @param sim_nets A list of simulated multiplex networks.
#' @param network_statistics A list of network statistics to compute for each network.
#' @param return_data A logical value indicating whether to return the computed statistics.
#' @param descriptive_labels Optional labels for the x-axis of the plot.
#'
#' @return If return_data is TRUE, a list containing the computed statistics for the simulated and observed networks.
#'         If return_data is FALSE, a ggplot object representing the comparison plot.
#'
#' @export
multiplex_gof_baseline <- function(dep_net, sim_nets, return_data = FALSE, descriptive_labels=NULL, ...){
    dep_lab <- names(dep_net)
    observed_stats <- as.data.frame(descriptive_stats(dep_net))
    observed_stats_df <- data.frame("var" = rownames(observed_stats), "sim_stats" = observed_stats[,1])
    sim_stats <- observed_stats
    num_sim <- length(sim_nets[[1]])
    stats <- list()
    for (i in 1:num_sim) {
        net <- lapply(dep_lab, function(x) sim_nets[[x]][[i]])
        names(net) <- dep_lab
        stats[[i]] <- descriptive_stats(net)
    }
    stats <- do.call(rbind, stats)
    basic_stats <- tidyr::pivot_longer(data.frame(stats), everything(), names_to="var", values_to="sim_stats")
    basic_stats$var <- factor(basic_stats$var, levels = colnames(stats))
    if (is.null(descriptive_labels)) {
        descriptive_labels <- levels(basic_stats$var)
    }
    p3 <- ggplot(basic_stats,aes(x = var, y=sim_stats, fill=var),show.legend = FALSE) +
 stat_boxplot(geom ='errorbar') +
    geom_boxplot() +
    scale_fill_manual(values=c(rep("#219ebc", 3), rep("#5AB1BB", 3), rep("#fb8500", 3), rep("#ffb703", 3)))+
    geom_point(data=observed_stats_df, mapping = aes(var,sim_stats), shape = 21, colour = "black", fill = "white", size=3) +
    theme_bw() + 
    theme(
      legend.position="none",
      plot.title = element_text(size=11),
      axis.text.y = element_text(size=15, color = "black"),
      axis.text.x = element_text(size=15, color = "black")) +
    scale_x_discrete(labels= descriptive_labels) + 
    xlab("") +ylab("") + coord_flip() +
    ggtitle("Descriptive Statistics Observed vs Simulated basic multiplex statistics")

    if (return_data) {
        return(list(stat_sim = basic_stats, stat_obs = observed_stats))
    } else {
        return(p3)
    }

}

#Currently this only works for biplex (t = 2) networks
multiplex_gof_random <- function(dep_net, sim_nets, return_data, descriptive_labels=NULL, ...){
    # if (length(dep_net) != 2) {
    #     stop("This function only works for biplex networks")
    # }
    dep_lab <- names(dep_net)
    observed_stats <- as.data.frame(descriptive_stats_Sigma_t3(dep_net))
    observed_stats_df <- data.frame("var" = rownames(observed_stats), "sim_stats" = observed_stats[,1])
    sim_stats <- observed_stats
    # post processing the generated results
    num_sim <- length(sim_nets[[1]])
    stats <- list()
    for (i in 1:num_sim) {
        net <- lapply(dep_lab, function(x) sim_nets[[x]][[i]])
        names(net) <- dep_lab
        stats[[i]] <- descriptive_stats_Sigma_t3(net)
    }
    stats <- data.frame(do.call(rbind, stats))
    basic_stats <- tidyr::pivot_longer(data.frame(stats), everything(), names_to="var", values_to="sim_stats")
    basic_stats$var <- factor(basic_stats$var, colnames(stats))
    iswithin <- sapply(strsplit(colnames(stats), "_"), function(x) x[1] == x[4])
    colors <- c("#219ebc", "#ffb703")
    p4 <- ggplot(basic_stats,aes(x = var, y=sim_stats, fill=var),show.legend = FALSE) +
    stat_boxplot(geom ='errorbar') +
        geom_boxplot() +
        scale_fill_manual(values=colors[iswithin + 1])+
        geom_point(data=observed_stats_df, mapping = aes(var,sim_stats), shape = 21, colour = "black", fill = "white", size=3) +
        theme_bw() + 
        theme(
        legend.position="none",
        plot.title = element_text(size=11),
        axis.text.y = element_text(size=15, color = "black"),
        axis.text.x = element_text(size=15, color = "black")) +
        #scale_x_discrete(labels= descriptive_labels) + 
        xlab("") +ylab("") + coord_flip()
        ggtitle("Descriptive Statistics Observed vs Simulated basic statistics") 


    if (return_data) {
        return(list(stat_sim = basic_stats, stat_obs = observed_stats))
    } else {
        return(p4)
    }
}

simulated_network_checks_single_layer <- function(sim_nets, dep_net, layer_lab, network_statistics_func, cumulative=TRUE) {
    num_sim <- length(sim_nets[[layer_lab]])
    obs_net <- dep_net[[layer_lab]]
    stat_obs <- network_statistics_func(obs_net)
    col_labs = names(stat_obs)
    stat_sim <- sapply(sim_nets[[layer_lab]], function(x) network_statistics_func(x))
    stat_sim <- matrix(stat_sim, ncol=num_sim)
	dimnames(stat_sim)[[1]] <- col_labs
	dimnames(stat_sim)[[2]] <-	1:num_sim
    stat_obs <- matrix(stat_obs, nrow=1)
    colnames(stat_obs) <- col_labs
	return(list(stat_sim = t(stat_sim), stat_obs = stat_obs))

}




#' simulated_network_checks Function
#'
#' This function performs predictive checks on a fitted model using simulated networks.
#'
#' @param model_obj The fitted model object.
#' @param network_statistics The network statistics function to be used for the checks.
#' @param num_sim The number of simulated networks to generate.
#' @param layer_lab The labels of the network layers to be used. If NULL, all network layers are used.
#' @param cumulative Logical indicating whether to calculate cumulative statistics.
#' @param return_data Logical indicating whether to return the simulated network statistics.
#' @param center Logical indicating whether to center the statistics.
#' @param scale Logical indicating whether to scale the statistics.
#' @param violin Logical indicating whether to plot the statistics using violin plots.
#' @param key The key for the violin plots.
#' @param perc The percentile for the key.
#' @param position The position of the key.
#' @param fontsize The font size of the plot.
#' @param ... Additional arguments to be passed to the network statistics function.
#'
#' @return If return_data is TRUE, a list of simulated network statistics is returned.
#'         Otherwise, the function plots the simulated network statistics.
#'
#' @details This function checks the goodness-of-fit of a fitted model by comparing the observed network statistics
#'          with the statistics calculated from the simulated networks. It supports both multiplex statistics and
#'          single layer statistics.
#'
#' @export
simulated_network_checks <- function(sim_nets, dep_net, network_statistics, layer_lab=NULL, cumulative=TRUE, return_data = FALSE, center=FALSE, scale=FALSE, violin=TRUE, key=NULL, perc=.05, position=4, fontsize=12, ...){
    print("Any missing value in the observed networks will be replaced with 0.")
    dep_net <- lapply(dep_net, function(x) {x[is.na(x)] <- 0; return(x)}) #remove NA values
    res_plots <- list()
    network_statistics_func <- match.fun(network_statistics)
    multiplex_gof_stats <- c("multiplex_gof_baseline", "multiplex_gof_random")

    if (is.null(layer_lab)) { #use all network layers
        layer_lab <- names(dep_net)
    }

    if (network_statistics %in% multiplex_gof_stats) { #workflows for multiplex statistics
        if (!is.null(layer_lab)) {
            print("all layers are used for multiplex statistics")
        }
        network_statistics_func(dep_net, sim_nets, return_data=return_data, ...)

    } else { #workflows for single layer statistics
        sim_stat_res <- vector("list", length(layer_lab))
        names(sim_stat_res) <- layer_lab

        for (layer in layer_lab) {
            sim_stat_res[[layer]] <- simulated_network_checks_single_layer(sim_nets, dep_net, layer, network_statistics_func, cumulative)
        }
        attr(sim_stat_res,"network_statistic_name") <- network_statistics
        attr(sim_stat_res,"key") <- names(sim_stat_res[[1]]$stat_obs) #maybe this is correct?
        
        if (return_data) {
            return(sim_stat_res)
        } else {
            res_plots = list()
            for (net_lab in layer_lab) {
                p <- gof_plot(sim_stat_res, net_lab, center=center, scale=scale, violin=violin, key=key, perc=perc, position=position, fontsize=fontsize, ...)
                res_plots[[net_lab]] <- p
            }
            gridExtra::grid.arrange(grobs = res_plots, ncol = length(layer_lab))
            return(res_plots)
        }
    }
    
}

#' Plot Goodness of Fit
#'
#' This function plots the goodness of fit for a Mp2 model
#'
#' @param sim_stat_res A list of simulated network statistics. Each element of the list should be a list with two elements: `stat_sim` and `stat_obs`. `stat_sim` is a matrix of simulated network statistics, and `stat_obs` is a matrix of observed network statistics.
#' @param net_lab The label of the network.
#' @param center Logical indicating whether to center the statistics.
#' @param scale Logical indicating whether to scale the statistics.
#' @param violin Logical indicating whether to plot violin plots.
#' @param key Optional vector of labels for the x-axis.
#' @param perc The percentage of simulations to use for the confidence intervals.
#' @param position The position of the observation labels.
#' @param fontsize The font size for the plot.
#' @param ... Additional arguments to be passed to the plotting functions.

#' @return None.
#' @import lattice
#' @export
gof_plot <- function (sim_stat_res, net_lab, center=FALSE, scale=FALSE, violin=TRUE, key=NULL, perc=.05, position=4, fontsize=12, ...)
{
	## require(lattice)
    args <- list(...)
	if (is.null(args$main)){
		main=paste("Goodness of Fit of",
				 attr(sim_stat_res,"network_statistic_name"), " for ", net_lab, sep=" ")
	} else {
		main=args$main
	}

	x <- sim_stat_res[[net_lab]]

	sims <- x$stat_sim
	obs <- x$stat_obs
    #colnames(obs) <- colnames(sims)
	itns <- nrow(sims)
	n.obs <- nrow(obs)

	screen <- sapply(1:ncol(obs),function(i){
						(sum(is.nan(rbind(sims,obs)[,i])) == 0) }) &
				(diag(var(rbind(sims,obs)))!=0)

	if (any((diag(var(rbind(sims,obs)))==0)))
	{	cat("Note: some statistics are not plotted because their variance is 0.\n")
		cat("This holds for the statistic")
		if (sum(diag(var(rbind(sims,obs)))==0) > 1){cat("s")}
		cat(": ")
		cat(paste(attr(x,"key")[which(diag(var(rbind(sims,obs)))==0)], sep=", "))
		cat(".\n")
	}

	sims <- sims[,screen, drop=FALSE]
	obsLabels <- round(obs[,screen, drop=FALSE],3)
    obs <- obs[,screen, drop=FALSE]

	sims.min <- apply(sims, 2, min)
	sims.max <- apply(sims, 2, max)
	sims.min <- pmin(sims.min, obs)
	sims.max <- pmax(sims.max, obs)

	if (center) {
		sims.median <- apply(sims, 2, median)
		sims <- sapply(1:ncol(sims), function(i)
					(sims[,i] - sims.median[i]) )
		obs <- matrix(sapply(1:ncol(sims), function(i)
							(obs[,i] - sims.median[i])), nrow=n.obs )
		sims.min <- sims.min - sims.median
		sims.max <- sims.max - sims.median
	} 
    
    if (scale) {
		sims.range <- sims.max - sims.min + 1e-6
		sims <- sapply(1:ncol(sims), function(i) sims[,i]/(sims.range[i]))
		obs <- matrix(sapply(1:ncol(sims), function(i) obs[,i]/(sims.range[i]))
				, nrow=n.obs )
		sims.min <- sims.min/sims.range
		sims.max <- sims.max/sims.range
	}

	ymin <- 1.05*min(sims.min) - 0.05*max(sims.max)
	ymax <- -0.05*min(sims.min) + 1.05*max(sims.max)

	if (is.null(args$ylab)){
		ylabel = "Statistic"
		if (center && scale) {
			ylabel = "Statistic (centered and scaled)"
		} else if (scale) {
			ylabel = "Statistic (scaled)"
		} else if (center) {
			ylabel = "Statistic (center)"
		} else {
			ylabel = "Statistic"
		}
	} else {
		ylabel = args$ylab
	} 
    
    if (is.null(args$xlab)) {
		xlabel = paste( paste("p:",
						collapse = " "), collapse = "\n")
	} else {
		xlabel = args$xlab
        }

	if (is.null(args$cex)) {
		cexpar <- par("cex")
    } else {
		cexpar <- args$cex
	}

	if (is.null(args$cex.axis)) {
    cexaxispar <- par("cex.axis")
    } else {
        cexaxispar <- args$cex.axis
    }

    if (is.null(args$cex.main)) {
        cexmainpar <- par("cex.main")
    } else {
        cexmainpar <- args$cex.main
    }

    if (is.null(args$cex.lab)) {
        cexlabpar <- par("cex.lab")
    } else {
        cexlabpar <- args$cex.lab
    }

    if (is.null(args$cex.sub)) {
        cexsubpar <- par("cex.sub")
    } else {
        cexsubpar <- args$cex.sub
    }

	xAxis <- (1:sum(screen))

	if (is.null(key)) {
		if (is.null(attr(x, "key"))){
			key=xAxis
		} else {
			key <- attr(x,"key")[screen]
		}
	} else {
		key <- key[screen] ## added 1.1-244
		if (length(key) != ncol(obs)) {
			stop("Key length does not match the number of variates.")
		}
	}

	br <- trellis.par.get("box.rectangle")
	br$col <- 1
	trellis.par.set("box.rectangle", br)
	bu <- trellis.par.get("box.umbrella")
	bu$col <- 1
	trellis.par.set("box.umbrella", bu)
	plot.symbol <- trellis.par.get("plot.symbol")
	plot.symbol$col <- "black"
	plot.symbol$pch <- 4
	plot.symbol$cex <- cexpar  # default 1
	trellis.par.set("plot.symbol", plot.symbol)

#	plot.axis <- trellis.par.get("axis.text")
#	plot.axis$cex <- cexaxispar # default 1
	trellis.par.set("axis.text", list(cex=cexaxispar))

#	plot.xlab <- trellis.par.get("par.xlab.text")
#	plot.xlab$cex <- cexlabpar # default 1
	trellis.par.set("par.xlab.text", list(cex=cexlabpar))

#	plot.ylab <- trellis.par.get("par.ylab.text")
#	plot.ylab$cex <- cexlabpar # default 1
	trellis.par.set("par.ylab.text", list(cex=cexlabpar))

#	plot.main <- trellis.par.get("par.main.text")
#	plot.main$cex <- cexmainpar # default 1.2
	trellis.par.set("par.main.text", list(cex=cexmainpar))

#	plot.font <- trellis.par.get("fontsize")
#	plot.font$text <- fontsize
	trellis.par.set("fontsize", list(text=fontsize))

	panelFunction <- function(..., x=x, y=y, box.ratio){
		ind.lower <- max( round(itns * perc/2), 1)
		ind.upper <- round(itns * (1-perc/2))
		yperc.lower <- sapply(1:ncol(sims), function(i)
					sort(sims[,i])[ind.lower]  )
		yperc.upper <- sapply(1:ncol(sims), function(i)
					sort(sims[,i])[ind.upper]  )
		if (violin) {
			panel.violin(x, y, box.ratio=box.ratio, col = "transparent", ...)
		}
		panel.bwplot(x, y, box.ratio=.1, fill = "gray", ...)
		panel.xyplot(xAxis, yperc.lower, lty=3, col = "gray", lwd=3, type="l",
				...)
		panel.xyplot(xAxis, yperc.upper, lty=3, col = "gray", lwd=3, type="l",
				...)
		for(i in 1:nrow(obs))
		{
			panel.xyplot(xAxis, obs[i,],  col="red", type="l", lwd=1, ...)
			panel.xyplot(xAxis, obs[i,],  col="red", type="p", lwd=3, pch=19,
					...)
			panel.text(xAxis, obs[i,], labels=obsLabels[i,], pos=position)
		}
	}

	p <- bwplot(as.numeric(sims)~rep(xAxis, each=itns), horizontal=FALSE,
			panel = panelFunction, xlab=xlabel, ylab=ylabel, ylim=c(ymin,ymax),
			scales=list(x=list(labels=key), y=list(draw=FALSE)),
			main=main)
    
    return(p)
}



