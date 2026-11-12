# rsconnect::setAccountInfo(name = "<account>", token = "<token>", secret = "<secret>")
remotes::install_github("junruidi/ActCR", upgrade = "never")
# Refresh the bundled paper from Writing/ so the deployed copy is the latest draft
dir.create("Shiny_tool/paper", showWarnings = FALSE)
stopifnot(file.copy("Writing/20261109_Paper_4_FDA_Xin.pdf", "Shiny_tool/paper/paper.pdf", overwrite = TRUE))
files <- c("app.R", "helpers.R", "reference.rds", "paper/paper.pdf",
           file.path("logo", list.files("Shiny_tool/logo")), file.path("examples", list.files("Shiny_tool/examples")))
rsconnect::deployApp("Shiny_tool", appFiles = files, appName = "fpca-residual-tool", appTitle = "Beyond the 24-Hour Rhythm", forceUpdate = TRUE, lint = FALSE)
