install.packages(c("usethis", "gitcreds"))

usethis::use_git_config(user.name = "DEIN_NAME", user.email = "DEINE_MAIL")

usethis::create_project("/project/ag-baer/Semaglutide")

usethis::use_git()

usethis::create_github_token()
gitcreds::gitcreds_set()

usethis::use_github(private = TRUE)

dir.create("scripts", showWarnings = FALSE)
