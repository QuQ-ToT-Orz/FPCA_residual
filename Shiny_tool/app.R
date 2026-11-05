library(shiny)
library(ggplot2)
source("helpers.R")

reference <- readRDS("reference.rds")
addResourcePath("logo", "logo")
hours <- (1:1440 - 1) / 60
title <- "Beyond the 24-Hour Rhythm: Daily Activity Residuals Patterns and Mortality"
minute_columns <- paste0("MIN", 1:1440)
scale_note <- "Activity is shown on a compressed scale, so brief bursts of very high activity do not dominate."
line_person <- "This person's pattern"
line_score <- "Reference average adjusted by this person's score"
line_average <- "Reference group average"
profile_names <- c("Group", "Participants", "Share of population, %", "Age, years", "Female, %", "BMI, kg/m2", "Average level", "Rhythm strength", "Peak time", "Daily pattern score")
term_names <- c(Age = "age", Gender = "gender", Income = "income", MESOR = "average activity level", Amplitude = "rhythm strength", Acrophase = "peak activity time", FPC2 = "daily pattern score")
# Personal variables in the cohort's risk model (every term not derived from the activity data); only these are asked for
person_terms <- function(cohort) setdiff(reference[[cohort]]$mortality$terms, c("MESOR", "Amplitude", "Acrophase", "FPC2"))
# The form has fields for age, gender and income only; a model needing anything else must get a field first
stopifnot(all(unlist(lapply(names(reference), person_terms)) %in% c("Age", "Gender", "Income")))
income_labels <- c("Below 1.3 times", "1.3 to below 3.5 times", "3.5 times or above", "Not reported")

ui <- fluidPage(
  title = title,
  tags$style(HTML("
    body { background: #f8f3ea; color: #1f2a2e; font-family: Georgia, serif; }
    .card { background: #fffcf6; border: 1px solid #e2d9cb; border-radius: 12px; padding: 18px; margin-bottom: 16px; }
    .note { color: #5d676c; font-size: 13px; }
    .btn-primary { background: #9b4d2d; border-color: #9b4d2d; }
    .card table { font-size: 13px; }
    summary { display: list-item; cursor: pointer; font-weight: bold; }
    .logos img { height: 40px; margin: 0 16px 16px 0; }")),
  h2(title),
  p("This research tool compares one person's wearable activity data with a reference group of US adults aged 20 and older ",
    "from the National Health and Nutrition Examination Survey (NHANES). It shows how the person's daily activity pattern compares, ",
    "which typical activity group it is closest to, and their estimated risk of cardiovascular death. ",
    "The reference models are fixed: uploaded data do not change them."),
  tags$ol(
    tags$li("Choose the device that recorded the data."),
    tags$li("Upload the activity file."),
    tags$li("Enter the person's details and click Calculate.")),
  fluidRow(
    column(4, div(class = "card",
      selectInput("cohort", "Device and NHANES reference group", c(
        "Hip, activity counts (2003-2006)" = "count",
        "Wrist, MIMS units (2011-2014)" = "mims")),
      fileInput("activity_file", "Activity file (.csv or .rds)", accept = c(".csv", ".rds")),
      p(class = "note", "One row per day, ideally 7 days, with one column per minute from MIN1 (midnight) to MIN1440. ",
        "Use the unit that matches the chosen device. If the file covers several people, add a SEQN column with their IDs. ",
        "Leave minutes when the device was not worn empty; they count as no activity."),
      selectInput("participant", "Participant ID", choices = NULL),
      uiOutput("person_inputs"),
      actionButton("run", "Calculate", class = "btn-primary"))),
    column(8, div(class = "card",
      conditionalPanel("input.run == 0", p(class = "note", "Results appear here after you click Calculate.")),
      tabsetPanel(
        tabPanel("Daily pattern", uiOutput("fpc2_text"), plotOutput("fpc2_plot", height = "380px")),
        tabPanel("Activity group", uiOutput("cluster_text"), plotOutput("cluster_plot", height = "340px"), tableOutput("cluster_table")),
        tabPanel("Cardiovascular risk", uiOutput("mortality_text")))))),
  div(class = "card", tags$details(
    tags$summary("Technical details"),
    tags$ul(
      tags$li("Activity is asinh-transformed. A cosinor model gives the 24-hour rhythm: average level (MESOR), strength (amplitude) and peak time (acrophase)."),
      tags$li("The daily pattern score is the second functional principal component (FPC2) of participant-average cosinor residuals, in cohort SD units."),
      tags$li("The hazard ratio is from a survey-weighted Cox model that also adjusts for FPC1 and FPC3."),
      tags$li(sprintf("Activity groups are from k-means on FPC scores of average activity, with K chosen by the gap statistic. Average silhouette is %.2f (hip) and %.2f (wrist), below the 0.25 that indicates substantial structure, so groups are descriptive.",
                      reference$count$cluster$silhouette, reference$mims$cluster$silhouette)),
      tags$li("Risk is from a survey-weighted logistic model with terms chosen by cross-validated AUC. The interval does not include error in the score or the term selection step.")))),
  div(class = "logos",
      tags$img(src = "logo/eu_marie_curie.png", alt = "EU Marie Curie"),
      tags$img(src = "logo/logo%20labda%202023%20horizontal%20large.png", alt = "Labda")))

server <- function(input, output, session) {
  activity <- reactive(read_table(req(input$activity_file)$datapath))
  observeEvent(activity(), updateSelectInput(session, "participant", choices = unique(activity()$SEQN)))

  # Values already entered are kept when the device changes
  output$person_inputs <- renderUI({
    terms <- person_terms(input$cohort)
    levels <- reference[[input$cohort]]$levels
    tagList(
      p(class = "note", "Only the details used by this device's risk model are asked for."),
      if ("Age" %in% terms) numericInput("age", "Age, years", isolate(input$age) %||% 60),
      if ("Gender" %in% terms) selectInput("gender", "Gender", levels$gender, isolate(input$gender)),
      if ("Income" %in% terms) selectInput("income", "Family income relative to the poverty line", setNames(levels$income, income_labels), isolate(input$income)))
  })

  ref <- eventReactive(input$run, reference[[input$cohort]])
  result <- eventReactive(input$run, {
    person <- list(Age = input$age, Gender = input$gender, Income = input$income)[person_terms(input$cohort)]
    problem <- if (is.null(input$activity_file)) "Please upload an activity file first." 
    else if (!all(minute_columns %in% names(activity()))) "The file needs one column per minute, named MIN1 to MIN1440."
    else if (any(lengths(person) == 0) || anyNA(unlist(person))) "Please fill in all the details."
    if (!is.null(problem)) {
      showNotification(problem, type = "error")
      req(FALSE)
    }
    x <- participant_matrix(activity(), input$participant)
    score_participant(x, reference[[input$cohort]], person)
  })

  output$fpc2_text <- renderUI({
    r <- result()
    tagList(
      h3(sprintf("Daily pattern score: %+.2f (higher than %.0f%% of the reference group)", r$fpc2, r$fpc2_percentile)),
      p("The score describes how activity is spread across the day once the regular 24-hour rest-activity rhythm is set aside. ",
        "A score of 0 is the reference group average, and one unit is one standard deviation. Higher scores mean relatively more ",
        "activity in the early morning and mid-afternoon, and relatively less in the mid-morning and late evening."),
      p(with(as.list(ref()$fpc2$hr), sprintf(
        "In this reference group, a score one unit higher was associated with %.2f times the rate of cardiovascular death (95%% CI %.2f to %.2f), after accounting for age, gender, race/ethnicity, income and the 24-hour rhythm. This is an association and does not show cause and effect.", hr, lo, hi))),
      p(sprintf("24-hour rhythm, from %d days of data: average level %.2f, strength %.2f, peak activity at %s.", r$days, r$mesor, r$amplitude, clock(r$acrophase))),
      if (r$days != 7) p(class = "note", sprintf("The reference groups used 7 days of data, so this result, based on %d days, is less comparable.", r$days)))
  })

  output$fpc2_plot <- renderPlot({
    r <- result()
    curves <- rbind(data.frame(hour = hours, value = r$residual, curve = line_person),
                    data.frame(hour = hours, value = r$fpc2_curve, curve = line_score),
                    data.frame(hour = hours, value = ref()$fpc2$mu, curve = line_average))
    ggplot(curves, aes(hour, value, colour = curve, linewidth = curve)) +
      geom_hline(yintercept = 0, colour = "grey70") + geom_line() + hour_axis +
      scale_colour_manual(values = setNames(c("grey60", "#9b4d2d", "#214b57"), c(line_person, line_score, line_average))) +
      scale_linewidth_manual(values = setNames(c(0.3, 0.9, 0.7), c(line_person, line_score, line_average))) +
      theme_bw(base_size = 13) + theme(legend.position = "bottom") + guides(colour = guide_legend(nrow = 2)) +
      labs(x = NULL, y = "Difference from the 24-hour rhythm", colour = NULL, linewidth = NULL, caption = scale_note)
  })

  output$cluster_text <- renderUI({
    tagList(
      h3(sprintf("Closest activity group: Group %d of %d", result()$cluster, nrow(ref()$cluster$centers))),
      p("People in the reference group were sorted into groups with similar average activity across the day. This shows which group ",
        "this person's data sit closest to. The groups overlap and are best read as points along a continuum, not as distinct types of people."))
  })

  output$cluster_plot <- renderPlot({
    r <- result()
    curves <- ref()$cluster$curves
    ggplot(curves, aes(hour, value)) +
      geom_line(aes(group = cluster), colour = "#c7b8ab") +
      geom_line(data = curves[curves$cluster == r$cluster, ], colour = "#9b4d2d", linewidth = 1) +
      geom_line(data = data.frame(hour = hours, value = r$curve), colour = "#214b57", linewidth = 0.3) +
      hour_axis + theme_bw(base_size = 13) +
      labs(x = NULL, y = "Average activity", caption = paste0("Dark blue: this person. Red-brown: average of the closest group. Tan: averages of the other groups.\n", scale_note))
  })

  output$cluster_table <- renderTable(setNames(ref()$cluster$profile, profile_names), digits = 1, striped = TRUE, caption = "Group averages in the reference data, weighted to the survey population. Average level and rhythm strength use the compressed activity scale.")

  output$mortality_text <- renderUI({
    r <- result()
    m <- ref()$mortality
    tagList(
      h3(sprintf("Estimated risk of cardiovascular death: %.1f%% (95%% CI %.1f%% to %.1f%%)", 100 * r$risk[1], 100 * r$risk[2], 100 * r$risk[3])),
      p(sprintf("This is the estimated chance of cardiovascular death during the reference group's follow-up (median %.1f years) for someone with these details and activity data. It is higher than the risk of %.0f%% of the reference group, whose average risk is %.1f%%.", ref()$followup, r$risk_percentile, 100 * m$mean)),
      p(paste0("The estimate uses: ", paste(term_names[m$terms], collapse = ", "), ".")),
      p("The daily pattern score was kept in both the hip and wrist models because it improved how well they predicted cardiovascular death in cross-validation within each reference group. The models have not been tested in outside data."),
      p("Follow-up length differs between the hip and wrist groups, so their risks cannot be compared directly."),
      p(class = "note", "The interval shows uncertainty in the average risk for people with these characteristics, not a range for this individual. This is a research tool; results are not a diagnosis or medical advice."))
  })
}

shinyApp(ui, server)
