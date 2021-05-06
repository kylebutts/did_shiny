#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(tidyverse)
library(did)
library(fixest)
library(fastDummies)

source("https://raw.githubusercontent.com/kylebutts/templates/master/ggplot_theme/theme_kyle.R")

ui <- fluidPage(
    theme = bslib::bs_theme(bootswatch = "united"),
    fluidRow(
        column(8, offset = 2,
            h1("Data Generating Process"),
            plotOutput("did")
        )
    ),
    fluidRow(
        column(2, offset = 2,
            h2("Treatment Group 1"),
            sliderInput("te1", "Treatment Effect:", 2, min = -10, max = 10),
            sliderInput("te_m1", "Treatment Effect Slope:", value = 0, min = -0.5, max = 0.5, step = 0.05),
            sliderInput("g1", "Treatment Date:", value = 2004, min = 2001, max = 2019, sep = "")
        ),
        column(2,
            h2("Treatment Group 2"),
            sliderInput("te2", "Treatment Effect:", 1, min = -10, max = 10),
            sliderInput("te_m2", "Treatment Effect Slope:", value = 0, min = -0.5, max = 0.5, step = 0.05),
            sliderInput("g2", "Treatment Date:", value = 2012, min = 2001, max = 2019, sep = "")
        ),
        column(2, offset = 1,
               h3("Presets"), br(),
               fluidRow(actionButton("btn_homogeneous", "Homogeneous Effects")), br(),
               fluidRow(actionButton("btn_het_level", "Heterogeneity in Levels")), br(),
               fluidRow(actionButton("btn_het_level_homo_slopes", "Heterogeneity in Levels (w/ Slopes)")), br(),
               fluidRow(actionButton("btn_het_slopes", "Heterogeneity in Slopes")), br(),
               fluidRow(actionButton("btn_het", "Heterogeneity in Levels and Slopes"))
        )
    ),
    fluidRow(hr(style="margin-bottom: 40px;")),
    fluidRow(
        column(6, 
               h1("TWFE Event Study"),
               p("This runs a TWFE Event-Study regression with all lead and legs besides -1 and the earliest relative time period:"),
               withMathJax("$$Y_{it} = \\alpha_i + \\alpha_t + \\sum_{k = -K+1}^{-2} \\gamma_{k}^{lead} D_{it}^k + \\sum_{k = 0}^L \\gamma_{k}^{lag} D_{it}^k + \\varepsilon_{it}$$")
        ),
        column(6, 
               h1("C&S Event Study"),
               p("This runs the R `did' program to estimate event-study parameters"),
        )
    ),
    fluidRow(
        column(6, plotOutput("es_twfe")),
        column(6, plotOutput("es_cs"))
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
            expand_grid(year = 2000:2020) %>%
            # Year FE
            group_by(year) %>% mutate(year_fe = rnorm(length(year), 0, 1)) %>% ungroup() %>%
            mutate(
                treat = g > year,
                rel_year = if_else(group == "Control", as.numeric(NA), as.numeric(year - g)),
                error = rnorm(n(), 0, 1),
                # Level Effect
                te = (group == "Group 1") * input$te1 * (year >= input$g1) + (group == "Group 2") * input$te2 * (year >= input$g2),
                # dynamic Effect
                te_dynamic = (group == "Group 1") * (year >= input$g1) * input$te_m1 * (year - input$g1) + (group == "Group 2") * (year >= input$g2) * input$te_m2 * (year - input$g2),
                dep_var = unit_fe + year_fe + te + te_dynamic + error
            )
    })
    
    # requires rel_year, estimate, se, ci_lower, ci_upper named as that
    make_es_plot <- function(es, te_true, estimator) {
        max_y <- max(es$ci_upper)
        
        es <- es %>% mutate(group = "Estimated Effect")
        te_true <- te_true %>% mutate(estimate = te_true, group = "True Effect")
        es <- bind_rows(es, te_true)
        
        max_y <- max(es$estimate)
        
        ggplot() + 
            # 0 effect
            geom_hline(yintercept = 0, linetype = "dashed") + 
            geom_vline(xintercept = -0.5, linetype = "dashed") +
            # Confidence Intervals
            geom_linerange(data = es, mapping = aes(x = rel_year, ymin = ci_lower, ymax = ci_upper), color = "grey30") +
            # Estimates
            geom_point(data = es, mapping = aes(x = rel_year, y = estimate, color = group), size = 2) + 
            # Label
            geom_label(data = data.frame(x = -0.5 - 0.1, y = max_y + 0.25, label = "Treatment Starts ▶"), label.size=NA,
                       mapping = aes(x = x, y = y, label = label), size = 5.5, hjust = 1, fontface = 2, inherit.aes = FALSE) +
            scale_x_continuous(breaks = -5:5, minor_breaks = NULL) +
            scale_y_continuous(minor_breaks = NULL) +
            scale_color_manual(values = c("Estimated Effect" = "#013ef5", "True Effect" = "#eb3f25")) +
            labs(x = "Relative Time", y = "Estimate", color = NULL, title = estimator) + 
            theme_kyle(base_size = 24) +
            theme(legend.position = "bottom")
    }
    
    get_true_effects <- function(df) {
        df %>%
            # Keep only treated units
            filter(g > 0) %>%
            group_by(rel_year) %>%
            summarize(te_true = mean(te + te_dynamic)) %>%
            filter(rel_year >= -5 & rel_year <= 5)
    }
    
    # Change second group treatment year options
    observeEvent(input$g1, {
        updateSliderInput(inputId = "g2", min = input$g1)
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
    observeEvent(input$btn_het_level_homo_slopes, {
        updateSliderInput(inputId = "te1", value = 2)
        updateSliderInput(inputId = "te2", value = input$te1 * 2)
        updateSliderInput(inputId = "te_m1", value = 0.1)
        updateSliderInput(inputId = "te_m2", value = 0.1)
    })
    observeEvent(input$btn_het_slopes, {
        updateSliderInput(inputId = "te1", value = input$te1)
        updateSliderInput(inputId = "te2", value = input$te1)
        updateSliderInput(inputId = "te_m1", value = 0.05)
        updateSliderInput(inputId = "te_m2", value = 0.2)
    })
    observeEvent(input$btn_het, {
        updateSliderInput(inputId = "te1", value = 2)
        updateSliderInput(inputId = "te2", value = input$te1 * 2)
        updateSliderInput(inputId = "te_m1", value = 0.05)
        updateSliderInput(inputId = "te_m2", value = 0.2)
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
    
    output$es_twfe <- renderPlot({
        df <- df()

        te_true <- get_true_effects(df)
        
        min_rel_year <- min(df$rel_year, na.rm=TRUE)
        
        df <- df %>% mutate(rel_year = if_else(is.na(rel_year), 1000, rel_year))
        
        formula <- as.formula(glue::glue("dep_var ~ i(rel_year, drop=c(-1, {min_rel_year}, 1000)) | unit + year"))
        
        mod <- fixest::feols(formula, data = df)

        es <- broom::tidy(mod) %>%
            filter(str_detect(term, "rel_year::")) %>% 
            rename(rel_year = term, se = std.error) %>%
            mutate(
                rel_year = as.numeric(str_remove(rel_year, "rel_year::")),
                ci_lower = estimate - 1.96 * se,
                ci_upper = estimate + 1.96 * se
            ) %>%
            filter(rel_year <= 5 & rel_year >= -5) %>%
            bind_rows(tibble(rel_year = -1, estimate = 0, se = 0, ci_lower = 0, ci_upper = 0))
        
        make_es_plot(es, te_true, "TWFE Event Study")
    })

    output$es_cs <- renderPlot({
        df <- df()

        te_true <- get_true_effects(df)
        
        mod <- did::att_gt(
                yname = "dep_var",
                tname = "year",
                idname = "unit",
                gname = "g",
                data = df
            ) %>%
            did::aggte(type = "dynamic")
        
        es <- broom::tidy(mod) %>%
            select(
                rel_year = event.time, estimate, se = std.error, ci_lower = point.conf.low, ci_upper = point.conf.high
            ) %>%
            filter(rel_year <= 5 & rel_year >= -5)
        
        make_es_plot(es, te_true, "Callaway and Sant'Anna")
    })
}

# Run the application 
# input <- list(g1 = 2004L, g2 = 2012L, te1 = 2.0, te2 = 2.0, te_m1=0.0, te_m2=0.0)

shinyApp(ui = ui, server = server)


