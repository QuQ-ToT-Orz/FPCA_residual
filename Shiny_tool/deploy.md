# 1. Rebuild reference.rds and examples/ so the Race edits take effect
Rscript Shiny_tool/build_reference.R

# 2. Optional: check the app locally first
Rscript -e 'shiny::runApp("Shiny_tool")'

# 3. Check that your shinyapps.io account is already set up on this machine
Rscript -e 'rsconnect::accounts()'

# 4. Deploy
Rscript Shiny_tool/deploy_shinyapps.R
