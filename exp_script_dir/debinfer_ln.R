options(echo=TRUE)
args <- commandArgs(trailingOnly = TRUE)

# args [1] - number of iterations
# args [2] - name of file with input data

library(deBInfer)

iter <- as.numeric (args [1])

# Load precompiled C code
dyn.load(paste("updated_equation", .Platform$dynlib.ext, sep = ""))


obs = read.table(args [2], header=TRUE)

# Reorder to put x7 to the correct column (in the external form of the equations it is the second variable, but internally it is the last one)
#indx7 <- which (colnames (obs) == "7")
#if (length (indx7) > 0) {
#    obs <- cbind (obs [, -indx7[1]], "x7"=obs[, "x7"])
#}


vars_presented = c()
for (i in 2:ncol(obs)) {
	vars_presented = c(vars_presented, as.numeric(gsub("x", "", colnames(obs)[i])))
}
print(vars_presented)
obs_model <- function(data, sim.data, samp) {
	res = 0
	for (i in 1:length(vars_presented)) {
		res = res + sum(dnorm(log(obs[,i+1]), mean = log(sim.data[,vars_presented[i] + 2]), sd = samp[[sprintf("sd.x%d", vars_presented[i])]], log = TRUE))
	}
	return(res)
}

userpr = read.table("parsed_userpriors.txt")
args = list()
# Add ODE parameters

non_fixed_params_qty = 0
args[[sprintf("a%d",0)]] = debinfer_par(name = sprintf("a%d",0), var.type = "de", fixed = TRUE, value = 0) #a0
for (i in 1:nrow(userpr)) {
	distr_type = userpr[i, 2]
	mean = log10(userpr[i, 3])
	sd = userpr[i, 4] / 3
	#mean = userpr[i, 3]
	#sd = userpr[i, 4]
	
	if (distr_type == "fixed") {
		args[[sprintf("a%d",i)]] = debinfer_par(name = sprintf("a%d",i), var.type = "de", fixed = TRUE, value = mean)
	}
	else {
		non_fixed_params_qty = non_fixed_params_qty + 1
		args[[sprintf("a%d",i)]] = debinfer_par(name = sprintf("a%d",i), var.type = "de", fixed = FALSE,
						  	  value = mean, prior = "norm", hypers = list(mean = mean, sd = sd),
							  prop.var = sd/10, samp.type="rw")
	}
}

# Variences for observed data
args[[sprintf("sd.x%d",0)]] = debinfer_par(name = sprintf("sd.x%d",0), var.type = "obs", fixed = TRUE, value = 0)  #x0
init_vals = read.table("initial_values_deBInfer_formatted.txt")
for (i in 1:nrow(init_vals)) {
	args[[sprintf("sd.x%d",i)]] = debinfer_par(name = sprintf("sd.x%d",i), var.type = "obs", fixed = FALSE,
						   value = 0.05, prior = "lnorm", hypers = list(meanlog = 0, sdlog = 1),
						   prop.var = c(3,4), samp.type = "rw-unif")
}

# Initial values
#init_vals = read.table("initial_values_deBInfer_formatted.txt")
args[[sprintf("x%d",0)]] <- debinfer_par(name = sprintf("x%d",0), var.type = "init", fixed = TRUE, value = 0) #x0
for (i in 1:nrow(init_vals)) {
	args[[sprintf("x%d",i)]] <- debinfer_par(name = sprintf("x%d",i), var.type = "init", fixed = TRUE, 
						 value = as.numeric(gsub(".*=", "", init_vals[i,1])))
}

mcmc.pars <- do.call(setup_debinfer, args)

mcmc_samples <- de_mcmc(N = iter, data = obs, de.model = "derivs", method="bdf",
			 jacfunc = "jac", dllname = "updated_equation",
			 initfunc = "initmod",
			 obs.model = obs_model, all.params = mcmc.pars,
			 Tmax = max(obs$time), data.times = obs$time, cnt = 1,
			 plot = FALSE, solver = "ode", verbose.mcmc = FALSE, verbose = FALSE)

samples <- as.data.frame (mcmc_samples$samples[,1:non_fixed_params_qty])
numlines <- dim (samples) [1]
samples [, "StepNum"] <- c (1:numlines)

cat ("#", file="deBInfer_output.txt", append=FALSE)
write.table(samples, "deBInfer_output.txt", sep="\t", row.names=FALSE, quote=FALSE, append=TRUE)


