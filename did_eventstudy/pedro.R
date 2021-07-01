# Load libraries and set baseline parameters
library(tidyverse)
library(lfe)
library(fastDummies)
library(fixest)
library(ggthemes)
library(did)
theme_set(theme_clean() + theme(plot.background = element_blank()))


## Set Parameters --------------------------------------------------------------
iseed  = 20201221
set.seed(iseed)

nrep <- 100  
true_mu <- 1


# All Groups Treated -----------------------------------------------------------

## Generate Data ---------------------------------------------------------------

# treated cohorts consist of 250 obs each, with the treatment effect still = true_mu on average
make_data <- function(nobs = 1000, 
					  nstates = 40) {
	
	# unit fixed effects (unobservd heterogeneity)
	unit <- tibble(
		unit = 1:nobs,
		# generate state
		state = sample(1:nstates, nobs, replace = TRUE),
		unit_fe = rnorm(nobs, state/5, 1),
		# generate instantaneous treatment effect
		#mu = rnorm(nobs, true_mu, 0.2)
		mu = true_mu
	)
	
	# year fixed effects (first part)
	year <- tibble(
		year = 1980:2010,
		year_fe = rnorm(length(year), 0, 1)
	)
	
	# Put the states into treatment groups
	treat_taus <- tibble(
		# sample the states randomly
		state = sample(1:nstates, nstates, replace = FALSE),
		# place the randomly sampled states into four treatment groups G_g
		cohort_year = sort(rep(c(1986, 1992, 1998, 2004), 10))
	)
	
	# make main dataset
	# full interaction of unit X year 
	expand_grid(unit = 1:nobs, year = 1980:2010) %>% 
		left_join(., unit) %>% 
		left_join(., year) %>% 
		left_join(., treat_taus) %>% 
		# make error term and get treatment indicators and treatment effects
		# Also get cohort specific trends (modify time FE)
		mutate(error = rnorm(nobs*31, 0, 1),
			   treat = ifelse(year >= cohort_year, 1, 0),
			   tau = ifelse(treat == 1, mu, 0),
			   year_fe = year_fe + 0.1*(year - cohort_year)
		) %>% 
		# calculate cumulative treatment effects
		group_by(unit) %>% 
		mutate(tau_cum = cumsum(tau)) %>% 
		ungroup() %>% 
		# calculate the dep variable
		mutate(dep_var = (2010 - cohort_year) + unit_fe + year_fe + tau_cum + error)
	
}



## All Groups Treated - DGP Visualization --------------------------------------

data <- make_data()

# plot
plot1 <- data %>% 
	ggplot(aes(x = year, y = dep_var, group = unit)) + 
	geom_line(alpha = 1/8, color = "grey") + 
	geom_line(data = data %>% 
			  	group_by(cohort_year, year) %>% 
			  	summarize(dep_var = mean(dep_var)),
			  aes(x = year, y = dep_var, group = factor(cohort_year),
			  	color = factor(cohort_year)),
			  size = 2) + 
	labs(x = "", y = "Value", color = "Treatment group   ") + 
	geom_vline(xintercept = 1986, color = '#E41A1C', size = 2) + 
	geom_vline(xintercept = 1992, color = '#377EB8', size = 2) + 
	geom_vline(xintercept = 1998, color = '#4DAF4A', size = 2) + 
	geom_vline(xintercept = 2004, color = '#984EA3', size = 2) + 
	scale_color_brewer(palette = 'Set1') + 
	theme(legend.position = 'bottom',
		  #legend.title = element_blank(), 
		  axis.title = element_text(size = 14),
		  axis.text = element_text(size = 12))  +
	ggtitle("One draw of the DGP with homogeneous effects across cohorts \n and with all groups being eventually treated")+
	theme(plot.title = element_text(hjust = 0.5, size=12))

plot1




## Binned TWFE ES & All Groups Treated -----------------------------------------

run_ES_DiD <- function(...) {
	# resimulate the data
	
	data <- make_data()
	
	data <- data %>% 
		mutate(
			rel_year = year - cohort_year,
			rel_year = case_when(
				rel_year < -5 ~ -6,
				rel_year > 5 ~ 6,
				TRUE ~ rel_year
			)
		)
	
	# estimate the model
	mod <- fixest::feols(dep_var ~ i(rel_year, drop=-1) | unit + year, data = data, cluster = ~state)
	
	# grab the estimates
	es <- broom::tidy(mod) %>%
		filter(str_detect(term, "rel_year::")) %>% 
		mutate(
			rel_year = as.numeric(str_remove(term, "rel_year::")),
		) %>%
		filter(rel_year <= 5 & rel_year >= -5) %>%
		select(t = rel_year, estimate)
	
	es
}

data_classical <- furrr::future_map_dfr(1:nrep, run_ES_DiD)

colors <- c("True Effect" = "red", "Estimated Effect" = "blue")

ES_plot_classical <- data_classical %>% 
	group_by(t) %>% 
	summarize(avg = mean(estimate),
			  sd = sd(estimate),
			  lower.ci = avg - 1.96*sd,
			  upper.ci = avg + 1.96*sd) %>% 
	bind_rows(tibble(t = -1, avg = 0, sd = 0, lower.ci = 0, upper.ci = 0)) %>% 
	mutate(true_tau = ifelse(t >= 0, (t + 1)* true_mu, 0)) %>% 
	ggplot(aes(x = t, y = avg)) + 
	#geom_linerange(aes(ymin = lower.ci, ymax = upper.ci), color = 'darkgrey', size = 2) + 
	geom_ribbon(aes(ymin = lower.ci, ymax = upper.ci), color = "lightgrey", alpha = 0.2) +
	geom_point(color = 'blue', size = 3) + 
	geom_line(aes(color = 'Estimated Effect'), size = 1) + 
	geom_line(aes(x = t, y = true_tau, color = 'True Effect'), linetype = "dashed", size = 2) + 
	geom_hline(yintercept = 0, linetype = "dashed") + 
	scale_x_continuous(breaks = -5:5) + 
	labs(x = "Relative Time", y = "Estimate") + 
	theme(axis.title = element_text(size = 14),
		  axis.text = element_text(size = 12)) +
	ggtitle("TWFE event-study regression with binned end-points")+
	scale_color_manual(values = colors) + 
	theme(plot.title = element_text(hjust = 0.5, size=12),
		  legend.position = "bottom", 
		  legend.title = element_blank())

ES_plot_classical



## Saturated TWFE ES & All Groups Treated --------------------------------------

run_ES_DiD_sat <- function(...) {
	
	
	# resimulate the data
	data <- make_data()
	
	# make dummy columns
	data <- data %>% 
		# make relative year indicator
		mutate(rel_year = year - cohort_year)
	
	min_year <- min(data$rel_year)
	
	formula <- as.formula(glue::glue("dep_var ~ i(rel_year, drop=c(-1, {min_year})) | unit + year"))
	
	mod <- fixest::feols(formula, data = data, cluster = ~state)
	
	es <- broom::tidy(mod) %>%
		filter(str_detect(term, "rel_year::")) %>% 
		mutate(
			rel_year = as.numeric(str_remove(term, "rel_year::")),
		) %>%
		filter(rel_year <= 5 & rel_year >= -5) %>%
		select(t = rel_year, estimate)
	
	es
}

data_sat <- map_dfr(1:nrep, run_ES_DiD_sat)

ES_plot_sat <- data_sat %>% 
	group_by(t) %>% 
	summarize(avg = mean(estimate),
			  sd = sd(estimate),
			  lower.ci = avg - 1.96*sd,
			  upper.ci = avg + 1.96*sd) %>% 
	bind_rows(tibble(t = -1, avg = 0, sd = 0, lower.ci = 0, upper.ci = 0)) %>% 
	mutate(true_tau = ifelse(t >= 0, (t + 1)* true_mu, 0)) %>% 
	ggplot(aes(x = t, y = avg)) + 
	#geom_linerange(aes(ymin = lower.ci, ymax = upper.ci), color = 'darkgrey', size = 2) + 
	geom_ribbon(aes(ymin = lower.ci, ymax = upper.ci), color = "lightgrey", alpha = 0.2) +
	geom_point(color = 'blue', size = 3) + 
	geom_line(aes(color = 'Estimated Effect'), size = 1) + 
	geom_line(aes(x = t, y = true_tau, color = 'True Effect'), linetype = "dashed", size = 2) + 
	geom_hline(yintercept = 0, linetype = "dashed") + 
	scale_x_continuous(breaks = -5:5) + 
	labs(x = "Relative Time", y = "Estimate") + 
	theme(axis.title = element_text(size = 14),
		  axis.text = element_text(size = 12)) +
	ggtitle("TWFE event-study regression with 'all' leads and lags")+
	scale_color_manual(values = colors) + 
	theme(plot.title = element_text(hjust = 0.5, size=12),
		  legend.position = "bottom", 
		  legend.title = element_blank())


ES_plot_sat


## C&S & All Groups Treated ----------------------------------------------------

# function to run ES DID
run_CS <- function(...) {
	
	# resimulate the data
	data <- make_data()
	
	mod <- did::att_gt(yname = "dep_var", 
					   tname = "year",
					   idname = "unit",
					   gname = "cohort_year",
					   control_group= "notyettreated",
					   bstrap = FALSE,
					   data = data,
					   print_details = FALSE)
	event_std <- did::aggte(mod, type = "dynamic")
	
	att.egt <- event_std$att.egt
	names(att.egt) <- event_std$egt
	
	# grab the obs we need
	broom::tidy(att.egt) %>% 
		filter(names %in% -5:5) %>% 
		mutate(t = -5:5, estimate = x) %>% 
		select(t, estimate)
}

data_CS <- map_dfr(1:nrep, run_CS)

ES_plot_CS <- data_CS %>% 
	group_by(t) %>% 
	summarize(avg = mean(estimate),
			  sd = sd(estimate),
			  lower.ci = avg - 1.96*sd,
			  upper.ci = avg + 1.96*sd) %>% 
	mutate(true_tau = ifelse(t >= 0, (t + 1)* true_mu, 0)) %>% 
	ggplot(aes(x = t, y = avg)) + 
	#geom_linerange(aes(ymin = lower.ci, ymax = upper.ci), color = 'darkgrey', size = 2) + 
	geom_ribbon(aes(ymin = lower.ci, ymax = upper.ci), color = "lightgrey", alpha = 0.2) +
	geom_point(color = 'blue', size = 3) + 
	geom_line(aes(color = 'Estimated Effect'), size = 1) + 
	geom_line(aes(x = t, y = true_tau, color = 'True Effect'), linetype = "dashed", size = 2) + 
	geom_hline(yintercept = 0, linetype = "dashed") + 
	scale_x_continuous(breaks = -5:5) + 
	labs(x = "Relative Time", y = "Estimate") + 
	theme(axis.title = element_text(size = 14),
		  axis.text = element_text(size = 12)) +
	ggtitle("Event-study-parameters estimated using Callaway and Sant'Anna (2020)\nComparison group: Not-yet-treated")+
	scale_color_manual(values = colors) + 
	theme(plot.title = element_text(hjust = 0.5, size=12),
		  legend.position = "bottom", 
		  legend.title = element_blank())

ES_plot_CS



# Never Treated Group ----------------------------------------------------------


## Generate Data ---------------------------------------------------------------

# treated cohorts consist of 250 obs each, with the treatment effect still = true_mu on average
make_data2 <- function(nobs = 1000, 
					   nstates = 40) {
	
	# unit fixed effects (unobservd heterogeneity)
	unit <- tibble(
		unit = 1:nobs,
		# generate state
		state = sample(1:nstates, nobs, replace = TRUE),
		unit_fe = rnorm(nobs, state/5, 1),
		# generate instantaneous treatment effect
		#mu = rnorm(nobs, true_mu, 0.2)
		mu = true_mu
	)
	
	# year fixed effects (first part)
	year <- tibble(
		year = 1980:2010,
		year_fe = rnorm(length(year), 0, 1)
	)
	
	# Put the states into treatment groups
	treat_taus <- tibble(
		# sample the states randomly
		state = sample(1:nstates, nstates, replace = FALSE),
		# place the randomly sampled states into 1\{t \ge g \}G_g
		cohort_year = sort(rep(c(1986, 1992, 1998, 2004), 10))
	)
	
	# make main dataset
	# full interaction of unit X year 
	expand_grid(unit = 1:nobs, year = 1980:2010) %>% 
		left_join(., unit) %>% 
		left_join(., year) %>% 
		left_join(., treat_taus) %>% 
		# make error term and get treatment indicators and treatment effects
		# Also get cohort specific trends (modify time FE)
		mutate(error = rnorm(nobs*31, 0, 1),
			   treat = ifelse((year >= cohort_year)* (cohort_year != 2004), 1, 0),
			   tau = ifelse(treat == 1, mu, 0),
			   year_fe = year_fe + 0.1*(year - cohort_year)
		) %>% 
		# calculate cumulative treatment effects
		group_by(unit) %>% 
		mutate(tau_cum = cumsum(tau)) %>% 
		ungroup() %>% 
		# calculate the dep variable
		mutate(dep_var = (2010 - cohort_year) + unit_fe + year_fe + tau_cum + error) %>%
		# Relabel 2004 cohort as never-treated
		mutate(cohort_year = ifelse(cohort_year == 2004, Inf, cohort_year))
	
}

## Never Treated Group - DGP Visualization -------------------------------------

data <- make_data2()

# plot
plot2 <- data %>% 
	ggplot(aes(x = year, y = dep_var, group = unit)) + 
	geom_line(alpha = 1/8, color = "grey") + 
	geom_line(data = data %>% 
			  	group_by(cohort_year, year) %>% 
			  	summarize(dep_var = mean(dep_var)),
			  aes(x = year, y = dep_var, group = factor(cohort_year),
			  	color = factor(cohort_year)),
			  size = 2) + 
	labs(x = "", y = "Value",  color = "Treatment group   ") + 
	geom_vline(xintercept = 1986, color = '#E41A1C', size = 2) + 
	geom_vline(xintercept = 1992, color = '#377EB8', size = 2) + 
	geom_vline(xintercept = 1998, color = '#4DAF4A', size = 2) + 
	theme(legend.position = 'bottom',
		  #legend.title = element_blank(), 
		  axis.title = element_text(size = 14),
		  axis.text = element_text(size = 12)) +
	scale_color_manual(labels = c("1986", "1992", "1998", "Never-treated"),
					   values = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3"))+
	ggtitle("One draw of the DGP with homogeneous effects across cohorts \n and with a never-treated group")+
	theme(plot.title = element_text(hjust = 0.5, size=12))

plot2 


## Binned TWFE & Never Treated Group -------------------------------------------


run_ES_DiD_never <- function(...) {
	
	
	# resimulate the data
	data <- make_data2()
	
	data <- data %>% 
		mutate(
			rel_year = year - cohort_year,
			rel_year = case_when(
				rel_year == -Inf ~ -Inf,
				rel_year < -5 ~ -6,
				rel_year > 5 ~ 6,
				TRUE ~ rel_year
			)
		)
	
	# estimate the model
	mod <- fixest::feols(dep_var ~ i(rel_year, drop=c(-1, -Inf)) | unit + year, data = data, cluster = ~state)
	
	# grab the estimates
	es <- broom::tidy(mod) %>%
		filter(str_detect(term, "rel_year::")) %>% 
		mutate(
			rel_year = as.numeric(str_remove(term, "rel_year::")),
		) %>%
		filter(rel_year <= 5 & rel_year >= -5) %>%
		select(t = rel_year, estimate)
	
	es
}

data_classical_never <- map_dfr(1:nrep, run_ES_DiD_never)

ES_plot_classical_never <- data_classical_never %>% 
	group_by(t) %>% 
	summarize(avg = mean(estimate),
			  sd = sd(estimate),
			  lower.ci = avg - 1.96*sd,
			  upper.ci = avg + 1.96*sd) %>% 
	bind_rows(tibble(t = -1, avg = 0, sd = 0, lower.ci = 0, upper.ci = 0)) %>% 
	mutate(true_tau = ifelse(t >= 0, (t + 1)* true_mu, 0)) %>% 
	ggplot(aes(x = t, y = avg)) + 
	#geom_linerange(aes(ymin = lower.ci, ymax = upper.ci), color = 'darkgrey', size = 2) + 
	geom_ribbon(aes(ymin = lower.ci, ymax = upper.ci), color = "lightgrey", alpha = 0.2) +
	geom_point(color = 'blue', size = 3) + 
	geom_line(aes(color = 'Estimated Effect'), size = 1) + 
	geom_line(aes(x = t, y = true_tau, color = 'True Effect'), linetype = "dashed", size = 2) + 
	geom_hline(yintercept = 0, linetype = "dashed") + 
	scale_x_continuous(breaks = -5:5) + 
	labs(x = "Relative Time", y = "Estimate") + 
	theme(axis.title = element_text(size = 14),
		  axis.text = element_text(size = 12))+
	ggtitle("TWFE event-study regression with binned end-points")+
	scale_color_manual(values = colors) + 
	theme(plot.title = element_text(hjust = 0.5, size=12),
		  legend.position = "bottom", 
		  legend.title = element_blank())

ES_plot_classical_never



## Saturated TWFE ES & Never Treated Group -------------------------------------

# function to run ES DID

run_ES_DiD_sat_never <- function(...) {
	
	# resimulate the data
	data <- make_data2()
	
	# make dummy columns
	data <- data %>% 
		# make relative year indicator
		mutate(rel_year = year - cohort_year)
	
	# get the minimum relative year - we need this to reindex
	min_year <- min(data$rel_year * (data$rel_year != -Inf), na.rm = T)
	
	
	formula <- as.formula(glue::glue("dep_var ~ i(rel_year, drop=c(-1, {min_year}, -Inf)) | unit + year"))
	
	mod <- fixest::feols(formula, data = data, cluster = ~state)
	
	es <- broom::tidy(mod) %>%
		filter(str_detect(term, "rel_year::")) %>% 
		mutate(
			rel_year = as.numeric(str_remove(term, "rel_year::")),
		) %>%
		filter(rel_year <= 5 & rel_year >= -5) %>%
		select(t = rel_year, estimate)
	
	es
	
	
}


data_sat_never <- map_dfr(1:nrep, run_ES_DiD_sat_never)

ES_plot_sat_never <- data_sat_never %>% 
	group_by(t) %>% 
	summarize(avg = mean(estimate),
			  sd = sd(estimate),
			  lower.ci = avg - 1.96*sd,
			  upper.ci = avg + 1.96*sd) %>% 
	bind_rows(tibble(t = -1, avg = 0, sd = 0, lower.ci = 0, upper.ci = 0)) %>% 
	mutate(true_tau = ifelse(t >= 0, (t + 1)* true_mu, 0)) %>% 
	ggplot(aes(x = t, y = avg)) + 
	#geom_linerange(aes(ymin = lower.ci, ymax = upper.ci), color = 'darkgrey', size = 2) + 
	geom_ribbon(aes(ymin = lower.ci, ymax = upper.ci), color = "lightgrey", alpha = 0.2) +
	geom_point(color = 'blue', size = 3) + 
	geom_line(aes(color = 'Estimated Effect'), size = 1) + 
	geom_line(aes(x = t, y = true_tau, color = 'True Effect'), linetype = "dashed", size = 2) + 
	geom_hline(yintercept = 0, linetype = "dashed") + 
	scale_x_continuous(breaks = -5:5) + 
	labs(x = "Relative Time", y = "Estimate") + 
	theme(axis.title = element_text(size = 14),
		  axis.text = element_text(size = 12))+
	ggtitle("TWFE event-study regression with 'all' leads and lags")+
	scale_color_manual(values = colors) + 
	theme(plot.title = element_text(hjust = 0.5, size=12),
		  legend.position = "bottom", 
		  legend.title = element_blank())

ES_plot_sat_never


## C&S & Never Treated Group ---------------------------------------------------
run_CS_never <- function(...) {
	
	# resimulate the data
	data <- make_data2()
	data$cohort_year[data$cohort_year==Inf] <- 0
	
	mod <- did::att_gt(yname = "dep_var", 
					   tname = "year",
					   idname = "unit",
					   gname = "cohort_year",
					   control_group= "never_treated",
					   bstrap = FALSE,
					   data = data,
					   print_details = FALSE)
	event_std <- did::aggte(mod, type = "dynamic")
	
	att.egt <- event_std$att.egt
	names(att.egt) <- event_std$egt
	
	# grab the obs we need
	broom::tidy(att.egt) %>% 
		filter(names %in% -5:5) %>% 
		mutate(t = -5:5, estimate = x) %>% 
		select(t, estimate)
}

data_CS_never <- map_dfr(1:nrep, run_CS_never)

ES_plot_CS_never <- data_CS_never %>% 
	group_by(t) %>% 
	summarize(avg = mean(estimate),
			  sd = sd(estimate),
			  lower.ci = avg - 1.96*sd,
			  upper.ci = avg + 1.96*sd) %>% 
	mutate(true_tau = ifelse(t >= 0, (t + 1)* true_mu, 0)) %>% 
	ggplot(aes(x = t, y = avg)) + 
	#geom_linerange(aes(ymin = lower.ci, ymax = upper.ci), color = 'darkgrey', size = 2) + 
	geom_ribbon(aes(ymin = lower.ci, ymax = upper.ci), color = "lightgrey", alpha = 0.2) +
	geom_point(color = 'blue', size = 3) + 
	geom_line(aes(color = 'Estimated Effect'), size = 1) + 
	geom_line(aes(x = t, y = true_tau, color = 'True Effect'), linetype = "dashed", size = 2) + 
	geom_hline(yintercept = 0, linetype = "dashed") + 
	scale_x_continuous(breaks = -5:5) + 
	labs(x = "Relative Time", y = "Estimate") + 
	theme(axis.title = element_text(size = 14),
		  axis.text = element_text(size = 12))+
	ggtitle("Event-study-parameters estimated using Callaway and Sant'Anna (2020)\nComparison group: Never-treated units")+
	scale_color_manual(values = colors) + 
	theme(plot.title = element_text(hjust = 0.5, size=12),
		  legend.position = "bottom", 
		  legend.title = element_blank())

ES_plot_CS_never


## C&S Not Yet Treated & Never Treated Group -----------------------------------

# function to run ES DID
run_CS_ny <- function(...) {
	
	# resimulate the data
	data <- make_data2()
	data$cohort_year[data$cohort_year==Inf] <- 0
	
	mod <- did::att_gt(yname = "dep_var", 
					   tname = "year",
					   idname = "unit",
					   gname = "cohort_year",
					   control_group= "notyettreated",
					   bstrap = FALSE,
					   data = data,
					   print_details = FALSE)
	event_std <- did::aggte(mod, type = "dynamic")
	
	att.egt <- event_std$att.egt
	names(att.egt) <- event_std$egt
	
	# grab the obs we need
	broom::tidy(att.egt) %>% 
		filter(names %in% -5:5) %>% 
		mutate(t = -5:5, estimate = x) %>% 
		select(t, estimate)
}

data_CS_ny <- map_dfr(1:nrep, run_CS_ny)

ES_plot_CS_ny <- data_CS_ny %>% 
	group_by(t) %>% 
	summarize(avg = mean(estimate),
			  sd = sd(estimate),
			  lower.ci = avg - 1.96*sd,
			  upper.ci = avg + 1.96*sd) %>% 
	mutate(true_tau = ifelse(t >= 0, (t + 1)* true_mu, 0)) %>% 
	ggplot(aes(x = t, y = avg)) + 
	#geom_linerange(aes(ymin = lower.ci, ymax = upper.ci), color = 'darkgrey', size = 2) + 
	geom_ribbon(aes(ymin = lower.ci, ymax = upper.ci), color = "lightgrey", alpha = 0.2) +
	geom_point(color = 'blue', size = 3) + 
	geom_line(aes(color = 'Estimated Effect'), size = 1) + 
	geom_line(aes(x = t, y = true_tau, color = 'True Effect'), linetype = "dashed", size = 2) + 
	geom_hline(yintercept = 0, linetype = "dashed") + 
	scale_x_continuous(breaks = -5:5) + 
	labs(x = "Relative Time", y = "Estimate") + 
	theme(axis.title = element_text(size = 14),
		  axis.text = element_text(size = 12))+
	ggtitle("Event-study-parameters estimated using Callaway and Sant'Anna (2020)\nComparison group: Not-yet-treated units")+
	scale_color_manual(values = colors) + 
	theme(plot.title = element_text(hjust = 0.5, size=12),
		  legend.position = "bottom", 
		  legend.title = element_blank())

ES_plot_CS_ny






# Heterogeneous Treatment Effects

## Generate data ---------------------------------------------------------------

# treated cohorts consist of 250 obs each, with the treatment effect still = true_mu on average
make_data3 <- function(nobs = 1000, 
					   nstates = 40) {
	
	# unit fixed effects (unobservd heterogeneity)
	unit <- tibble(
		unit = 1:nobs,
		# generate state
		state = sample(1:nstates, nobs, replace = TRUE),
		unit_fe = rnorm(nobs, state/5, 1),
		# generate instantaneous treatment effect
		#mu = rnorm(nobs, true_mu, 0.2)
		mu = true_mu
	)
	
	# year fixed effects (first part)
	year <- tibble(
		year = 1980:2010,
		year_fe = rnorm(length(year), 0, 1)
	)
	
	# Put the states into treatment groups
	treat_taus <- tibble(
		# sample the states randomly
		state = sample(1:nstates, nstates, replace = FALSE),
		# place the randomly sampled states into 1\{t \ge g \}G_g
		cohort_year = sort(rep(c(1986, 1992, 1998, 2004), 10))
	)
	
	# make main dataset
	# full interaction of unit X year 
	expand_grid(unit = 1:nobs, year = 1980:2010) %>% 
		left_join(., unit) %>% 
		left_join(., year) %>% 
		left_join(., treat_taus) %>% 
		# make error term and get treatment indicators and treatment effects
		# Also get cohort specific trends (modify time FE)
		mutate(error = rnorm(nobs*31, 0, 1),
			   treat = ifelse((year >= cohort_year)* (cohort_year != 2004), 1, 0),
			   mu = ifelse(cohort_year==1992, 2, ifelse(cohort_year==1998, 1, 3)),
			   tau = ifelse(treat == 1, mu, 0),
			   year_fe = year_fe + 0.1*(year - cohort_year)
		) %>% 
		# calculate cumulative treatment effects
		group_by(unit) %>% 
		mutate(tau_cum = cumsum(tau)) %>% 
		ungroup() %>% 
		# calculate the dep variable
		mutate(dep_var = (2010 - cohort_year) + unit_fe + year_fe + tau_cum + error) %>%
		# Relabel 2004 cohort as never-treated
		mutate(cohort_year = ifelse(cohort_year == 2004, Inf, cohort_year))
	
}


## Heterogeneous -- Data Visualization -----------------------------------------

data <- make_data3()

# plot
plot3 <- data %>% 
	ggplot(aes(x = year, y = dep_var, group = unit)) + 
	geom_line(alpha = 1/8, color = "grey") + 
	geom_line(data = data %>% 
			  	group_by(cohort_year, year) %>% 
			  	summarize(dep_var = mean(dep_var)),
			  aes(x = year, y = dep_var, group = factor(cohort_year),
			  	color = factor(cohort_year)),
			  size = 2) + 
	labs(x = "", y = "Value",  color = "Treatment group   ") + 
	geom_vline(xintercept = 1986, color = '#E41A1C', size = 2) + 
	geom_vline(xintercept = 1992, color = '#377EB8', size = 2) + 
	geom_vline(xintercept = 1998, color = '#4DAF4A', size = 2) + 
	theme(legend.position = 'bottom',
		  #legend.title = element_blank(), 
		  axis.title = element_text(size = 14),
		  axis.text = element_text(size = 12)) +
	scale_color_manual(labels = c("1986", "1992", "1998", "Never-treated"),
					   values = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3")) +
	ggtitle("One draw of the DGP with heterogeneous treatment effect dynamics across cohorts \n and with a never-treated group")+
	theme(plot.title = element_text(hjust = 0.5, size=12))

plot3




## Binned TWFE ES - Heterogeneous Effects --------------------------------------

keepvars <- c("`rel_year_-5`",  "`rel_year_-4`",  "`rel_year_-3`",  "`rel_year_-2`",
			  "rel_year_0", "rel_year_1", "rel_year_2", "rel_year_3", "rel_year_4", "rel_year_5")


run_ES_DiD_never_het <- function(...) {
	
	
	# resimulate the data
	data <- make_data3()
	
	data <- data %>% 
		mutate(
			rel_year = year - cohort_year,
			rel_year = case_when(
				rel_year == -Inf ~ -Inf,
				rel_year < -5 ~ -6,
				rel_year > 5 ~ 6,
				TRUE ~ rel_year
			)
		)
	
	# estimate the model
	mod <- fixest::feols(dep_var ~ i(rel_year, drop=c(-1, -Inf)) | unit + year, data = data, cluster = ~state)
	
	# grab the estimates
	es <- broom::tidy(mod) %>%
		filter(str_detect(term, "rel_year::")) %>% 
		mutate(
			rel_year = as.numeric(str_remove(term, "rel_year::")),
		) %>%
		filter(rel_year <= 5 & rel_year >= -5) %>%
		select(t = rel_year, estimate)
	
	es
}


data_classical_never_het <- map_dfr(1:nrep, run_ES_DiD_never_het)

ES_plot_classical_never_het <- data_classical_never_het %>% 
	group_by(t) %>% 
	summarize(avg = mean(estimate),
			  sd = sd(estimate),
			  lower.ci = avg - 1.96*sd,
			  upper.ci = avg + 1.96*sd) %>% 
	bind_rows(tibble(t = -1, avg = 0, sd = 0, lower.ci = 0, upper.ci = 0)) %>% 
	mutate(true_tau = ifelse(t >= 0, (t + 1)* 2, 0)) %>% 
	ggplot(aes(x = t, y = avg)) + 
	#geom_linerange(aes(ymin = lower.ci, ymax = upper.ci), color = 'darkgrey', size = 2) + 
	geom_ribbon(aes(ymin = lower.ci, ymax = upper.ci), color = "lightgrey", alpha = 0.2) +
	geom_point(color = 'blue', size = 3) + 
	geom_line(aes(color = 'Estimated Effect'), size = 1) + 
	geom_line(aes(x = t, y = true_tau, color = 'True Effect'), linetype = "dashed", size = 2) + 
	geom_hline(yintercept = 0, linetype = "dashed") + 
	scale_x_continuous(breaks = -5:5) + 
	labs(x = "Relative Time", y = "Estimate") + 
	theme(axis.title = element_text(size = 14),
		  axis.text = element_text(size = 12))+
	ggtitle("TWFE event-study regression with binned end-points")+
	scale_color_manual(values = colors) + 
	theme(plot.title = element_text(hjust = 0.5, size=12),
		  legend.position = "bottom", 
		  legend.title = element_blank())

ES_plot_classical_never_het



## Saturated TWFE ES & Heterogeneous Effects -----------------------------------

# function to run ES DID
run_ES_DiD_sat_never_het <- function(...) {
	
	# resimulate the data
	data <- make_data3()
	
	# make dummy columns
	data <- data %>% 
		# make relative year indicator
		mutate(rel_year = year - cohort_year)
	
	# get the minimum relative year - we need this to reindex
	min_year <- min(data$rel_year * (data$rel_year != -Inf), na.rm = T)
	
	
	formula <- as.formula(glue::glue("dep_var ~ i(rel_year, drop=c(-1, {min_year}, -Inf)) | unit + year"))
	
	mod <- fixest::feols(formula, data = data, cluster = ~state)
	
	es <- broom::tidy(mod) %>%
		filter(str_detect(term, "rel_year::")) %>% 
		mutate(
			rel_year = as.numeric(str_remove(term, "rel_year::")),
		) %>%
		filter(rel_year <= 5 & rel_year >= -5) %>%
		select(t = rel_year, estimate)
	
	es
	
}

data_sat_never_het <- map_dfr(1:nrep, run_ES_DiD_sat_never_het)

ES_plot_sat_never_het <- data_sat_never_het %>% 
	group_by(t) %>% 
	summarize(avg = mean(estimate),
			  sd = sd(estimate),
			  lower.ci = avg - 1.96*sd,
			  upper.ci = avg + 1.96*sd) %>% 
	bind_rows(tibble(t = -1, avg = 0, sd = 0, lower.ci = 0, upper.ci = 0)) %>% 
	mutate(true_tau = ifelse(t >= 0, (t + 1)* 2, 0)) %>% 
	ggplot(aes(x = t, y = avg)) + 
	#geom_linerange(aes(ymin = lower.ci, ymax = upper.ci), color = 'darkgrey', size = 2) + 
	geom_ribbon(aes(ymin = lower.ci, ymax = upper.ci), color = "lightgrey", alpha = 0.2) +
	geom_point(color = 'blue', size = 3) + 
	geom_line(aes(color = 'Estimated Effect'), size = 1) + 
	geom_line(aes(x = t, y = true_tau, color = 'True Effect'), linetype = "dashed", size = 2) + 
	geom_hline(yintercept = 0, linetype = "dashed") + 
	scale_x_continuous(breaks = -5:5) + 
	labs(x = "Relative Time", y = "Estimate") + 
	theme(axis.title = element_text(size = 14),
		  axis.text = element_text(size = 12))+
	ggtitle("TWFE event-study regression with 'all' leads and lags")+
	scale_color_manual(values = colors) + 
	theme(plot.title = element_text(hjust = 0.5, size=12),
		  legend.position = "bottom", 
		  legend.title = element_blank())

ES_plot_sat_never_het



## C&S & Heterogeneous Effects -------------------------------------------------

# function to run ES DID
run_CS_never_het <- function(...) {
	
	# resimulate the data
	data <- make_data3()
	data$cohort_year[data$cohort_year==Inf] <- 0
	
	mod <- did::att_gt(yname = "dep_var", 
					   tname = "year",
					   idname = "unit",
					   gname = "cohort_year",
					   control_group= "never_treated",
					   bstrap = FALSE,
					   data = data,
					   print_details = FALSE)
	event_std <- did::aggte(mod, type = "dynamic")
	
	att.egt <- event_std$att.egt
	names(att.egt) <- event_std$egt
	
	# grab the obs we need
	broom::tidy(att.egt) %>% 
		filter(names %in% -5:5) %>% 
		mutate(t = -5:5, estimate = x) %>% 
		select(t, estimate)
}

data_CS_never_het <- map_dfr(1:nrep, run_CS_never_het)

ES_plot_CS_never_het <- data_CS_never_het %>% 
	group_by(t) %>% 
	summarize(avg = mean(estimate),
			  sd = sd(estimate),
			  lower.ci = avg - 1.96*sd,
			  upper.ci = avg + 1.96*sd) %>% 
	mutate(true_tau = ifelse(t >= 0, (t + 1)* 2, 0)) %>% 
	ggplot(aes(x = t, y = avg)) + 
	#geom_linerange(aes(ymin = lower.ci, ymax = upper.ci), color = 'darkgrey', size = 2) + 
	geom_ribbon(aes(ymin = lower.ci, ymax = upper.ci), color = "lightgrey", alpha = 0.2) +
	geom_point(color = 'blue', size = 3) + 
	geom_line(aes(color = 'Estimated Effect'), size = 1) + 
	geom_line(aes(x = t, y = true_tau, color = 'True Effect'), linetype = "dashed", size = 2) + 
	geom_hline(yintercept = 0, linetype = "dashed") + 
	scale_x_continuous(breaks = -5:5) + 
	labs(x = "Relative Time", y = "Estimate") + 
	theme(axis.title = element_text(size = 14),
		  axis.text = element_text(size = 12))+
	ggtitle("Event-study-parameters estimated using Callaway and Sant'Anna (2020)\nComparison group: Never-treated units")+
	scale_color_manual(values = colors) + 
	theme(plot.title = element_text(hjust = 0.5, size=12),
		  legend.position = "bottom", 
		  legend.title = element_blank())

ES_plot_CS_never_het



## C&S Not Yet Treated & Heterogeneous Effects ---------------------------------

# function to run ES DID
run_CS_ny_het <- function(...) {
	
	# resimulate the data
	data <- make_data3()
	data$cohort_year[data$cohort_year==Inf] <- 0
	
	mod <- did::att_gt(yname = "dep_var", 
					   tname = "year",
					   idname = "unit",
					   gname = "cohort_year",
					   control_group= "notyettreated",
					   bstrap = FALSE,
					   data = data,
					   print_details = FALSE)
	event_std <- did::aggte(mod, type = "dynamic")
	
	att.egt <- event_std$att.egt
	names(att.egt) <- event_std$egt
	
	# grab the obs we need
	broom::tidy(att.egt) %>% 
		filter(names %in% -5:5) %>% 
		mutate(t = -5:5, estimate = x) %>% 
		select(t, estimate)
}

data_CS_ny_het <- map_dfr(1:nrep, run_CS_ny_het)

ES_plot_CS_ny_het <- data_CS_ny_het %>% 
	group_by(t) %>% 
	summarize(avg = mean(estimate),
			  sd = sd(estimate),
			  lower.ci = avg - 1.96*sd,
			  upper.ci = avg + 1.96*sd) %>% 
	mutate(true_tau = ifelse(t >= 0, (t + 1)* 2, 0)) %>% 
	ggplot(aes(x = t, y = avg)) + 
	#geom_linerange(aes(ymin = lower.ci, ymax = upper.ci), color = 'darkgrey', size = 2) + 
	geom_ribbon(aes(ymin = lower.ci, ymax = upper.ci), color = "lightgrey", alpha = 0.2) +
	geom_point(color = 'blue', size = 3) + 
	geom_line(aes(color = 'Estimated Effect'), size = 1) + 
	geom_line(aes(x = t, y = true_tau, color = 'True Effect'), linetype = "dashed", size = 2) + 
	geom_hline(yintercept = 0, linetype = "dashed") + 
	scale_x_continuous(breaks = -5:5) + 
	labs(x = "Relative Time", y = "Estimate") + 
	theme(axis.title = element_text(size = 14),
		  axis.text = element_text(size = 12))+
	ggtitle("Event-study-parameters estimated using Callaway and Sant'Anna (2020)\nComparison group: Not-yet-treated units")+
	scale_color_manual(values = colors) + 
	theme(plot.title = element_text(hjust = 0.5, size=12),
		  legend.position = "bottom", 
		  legend.title = element_blank())

ES_plot_CS_ny_het







