# rsconnect::setAccountInfo(name = "<account>", token = "<token>", secret = "<secret>")
remotes::install_github("junruidi/ActCR", upgrade = "never")
files <- c("app.R", "helpers.R", "reference.rds",
           file.path("logo", list.files("Shiny_tool/logo")), file.path("examples", list.files("Shiny_tool/examples")))
rsconnect::deployApp("Shiny_tool", appFiles = files, appName = "fpca-residual-tool", appTitle = "Beyond the 24-Hour Rhythm", forceUpdate = TRUE, lint = FALSE)
