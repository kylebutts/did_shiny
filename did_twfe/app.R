## -----------------------------------------------------------------------------
## app.R
## Kyle Butts, CU Boulder Economics 
## 
## This is a shiny application to highlight the pitfals with TWFE DiD.
## -----------------------------------------------------------------------------

library(shiny)
library(tidyverse)
library(fixest)
library(did)

source("https://raw.githubusercontent.com/kylebutts/templates/master/ggplot_theme/theme_kyle.R")

ui <- fluidPage(
	theme = bslib::bs_theme(bootswatch = "united"),
	fluidRow(
		column(10, offset = 1,
			   h1("Data Generating Process"),
			   fluidRow(
					column(9, plotOutput("did")),
					column(3, tableOutput("did_estimates"))
			   )
		)
	),
	fluidRow(
		column(2, offset = 2,
			   h2("Treatment Group 1"),
			   sliderInput("te1", "Treatment Effect:", 2, min = -10, max = 10),
			   sliderInput("te_m1", "Treatment Effect Slope:", value = 0.3, min = -1, max = 1, step = 0.05),
			   sliderInput("g1", "Treatment Date:", value = 2004, min = 2001, max = 2019, sep = "")
		),
		column(2,
			   h2("Treatment Group 2"),
			   sliderInput("te2", "Treatment Effect:", 1, min = -10, max = 10),
			   sliderInput("te_m2", "Treatment Effect Slope:", value = 0.15, min = -1, max = 1, step = 0.05),
			   sliderInput("g2", "Treatment Date:", value = 2012, min = 2001, max = 2019, sep = "")
		),
		column(2,
			   h2("Panel"),
			   sliderInput("panel", "Panel Years:", c(2000, 2020), min = 1980, max = 2030, step = 5, sep = ""),
		),
		column(2,
			   h3("Presets"), br(),
			   fluidRow(actionButton("btn_homogeneous", "Homogeneous Effects")), br(),
			   fluidRow(actionButton("btn_het_levels", "Heterogeneity in Levels")), br(),
			   fluidRow(actionButton("btn_het_levels_homo_slopes", "Heterogeneity in Levels (w/ Slopes)")), br(),
			   fluidRow(actionButton("btn_het_slopes", "Heterogeneity in Slopes")), br(),
			   fluidRow(actionButton("btn_het", "Heterogeneity in Levels and Slopes"))
		)
	)
)

server <- function(input, output, session) {
	
	# Gen Data
	df <- reactive({
		tibble(unit = 1:1000) %>%
			mutate(
				state = sample(1:40, n(), replace = TRUE),
				unit_fe = rnorm(n(), state/5, 1),
				group = runif(n()),
				group = case_when(
					group < 0.33 ~ "Group 1",
					group < 0.66 ~ "Group 2",
					TRUE ~ "Control"
				),
				g = case_when(
					group == "Group 1" ~ input$g1,
					group == "Group 2" ~ input$g2,
					TRUE ~ 0L
				)
			) %>% 
			expand_grid(year = input$panel[1]:input$panel[2]) %>%
			# Year FE
			group_by(year) %>% mutate(year_fe = rnorm(length(year), 0, 1)) %>% ungroup() %>%
			mutate(
				treat = if_else(group == "Control", FALSE, year >= g),
				rel_year = if_else(group == "Control", Inf, as.numeric(year - g)),
				rel_year_binned = case_when(
					rel_year == Inf ~ Inf,
					rel_year <= -6 ~ -6,
					rel_year >= 6 ~ 6,
					TRUE ~ rel_year
				),
				error = rnorm(n(), 0, 1),
				# Level Effect
				te = (group == "Group 1") * input$te1 * (year >= input$g1) + (group == "Group 2") * input$te2 * (year >= input$g2),
				# dynamic Effect
				te_dynamic = (group == "Group 1") * (year >= input$g1) * input$te_m1 * (year - input$g1) + (group == "Group 2") * (year >= input$g2) * input$te_m2 * (year - input$g2),
				dep_var = unit_fe + year_fe + te + te_dynamic + error
			)
	})
	
	observeEvent(input$panel, {
		updateSliderInput(inputId = "g1", min = input$panel[1] + 2, max = input$panel[2] - 2)
		updateSliderInput(inputId = "g2", min = input$panel[1] + 2, max = input$panel[2] - 2)
	})
	
	# Presets
	observeEvent(input$btn_homogeneous, {
		updateSliderInput(inputId = "te1", value = input$te1)
		updateSliderInput(inputId = "te2", value = input$te1)
		updateSliderInput(inputId = "te_m1", value = 0)
		updateSliderInput(inputId = "te_m2", value = 0)
	})
	observeEvent(input$btn_het_levels, {
		updateSliderInput(inputId = "te1", value = 2)
		updateSliderInput(inputId = "te2", value = input$te1 * 2)
		updateSliderInput(inputId = "te_m1", value = 0)
		updateSliderInput(inputId = "te_m2", value = 0)
	})
	observeEvent(input$btn_het_levels_homo_slopes, {
		updateSliderInput(inputId = "te1", value = 2)
		updateSliderInput(inputId = "te2", value = input$te1 * 2)
		updateSliderInput(inputId = "te_m1", value = 0.15)
		updateSliderInput(inputId = "te_m2", value = 0.15)
	})
	observeEvent(input$btn_het_slopes, {
		updateSliderInput(inputId = "te1", value = input$te1)
		updateSliderInput(inputId = "te2", value = input$te1)
		updateSliderInput(inputId = "te_m1", value = 0.3)
		updateSliderInput(inputId = "te_m2", value = 0.15)
	})
	observeEvent(input$btn_het, {
		updateSliderInput(inputId = "te1", value = 2)
		updateSliderInput(inputId = "te2", value = input$te1 * 2)
		updateSliderInput(inputId = "te_m1", value = 0.3)
		updateSliderInput(inputId = "te_m2", value = 0.15)
	})
	
	
	output$did <- renderPlot({
		df <- df()
		
		df_avg <- df %>% group_by(group, year) %>% summarize(dep_var = mean(dep_var))
		
		max_y <- max(df_avg$dep_var)
		
		ggplot() + 
			geom_line(data = df_avg, mapping = aes(y = dep_var, x = year, color = group), size = 1.5) +
			geom_vline(xintercept = input$g1 - 0.5, linetype = "dashed") + 
			geom_vline(xintercept = input$g2 - 0.5, linetype = "dashed") +
			theme_kyle(base_size = 24) +
			theme(legend.position = "bottom") +
			labs(y = "Outcome", x = "Year", color = "Treatment Cohort") + 
			scale_y_continuous(expand = expansion(add = .5)) + 
			scale_color_manual(values = c("Control" = "#8e549f", "Group 1" = "#d2382c", "Group 2" = "#497eb3")) +
			geom_label(data = data.frame(x = input$g1 - 0.4, y = max_y + 0.5, label = "◀ Treatment Starts"), label.size=NA,
					   mapping = aes(x = x, y = y, label = label), size = 4.23, hjust = 0L, fontface = 2, inherit.aes = FALSE) +
			geom_label(data = data.frame(x = input$g2 - 0.4, y = max_y + 0.75, label = "◀ Treatment Starts"), label.size=NA,
					   mapping = aes(x = x, y = y, label = label), size = 4.23, hjust = 0L, fontface = 2, inherit.aes = FALSE)
	}, res = 96)
	
	output$did_estimates <- renderTable({
		df <- df()
		
		te_true <- df %>% filter(treat == TRUE) %>% mutate(effect = te + te_dynamic) %>% pull(effect) %>% mean()
			
		te_twfe <- fixest::feols(dep_var ~ treat | unit + year, data = df) %>% coefficients() %>% .[["treatTRUE"]] 
		
		te_did <- did::att_gt(
			yname = "dep_var", tname = "year", idname = "unit", gname = "g", data = df
		) %>% did::aggte(type = "simple") %>% broom::tidy() %>% pull(estimate)
		
		tribble(
			~Estimator, ~Effect,
			"True Treatment Effect", te_true,
			"TWFE Estimate", te_twfe,
			"C&S Estimate", te_did
		)
	})
}

# Run the application 
# For Testing
# input <- list(g1 = 2004L, g2 = 2012L, te1 = 2.0, te2 = 2.0, te_m1=1, te_m2=0.5, panel = c(2000, 2020))

shinyApp(ui = ui, server = server)


