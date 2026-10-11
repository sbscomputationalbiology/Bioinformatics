options(repos = c(CRAN = "https://packagemanager.posit.co/cran/__linux__/manylinux_2_28/latest"))
renv::install(c("usethis", "gitcreds"))

renv::snapshot(packages = c("usethis", "gitcreds"))

renv::status()

# Tell Git who you are (use the email on your GitHub account)
usethis::use_git_config(user.name = "sbscomputationalbiology", user.email = "sbscomputationalbiology@outlook.com")

# Create a personal access token (opens GitHub in your browser; accept the defaults and copy the token)
usethis::create_github_token()

# Store it (paste the token when prompted)
gitcreds::gitcreds_set()

# Check that everything is wired up
usethis::git_sitrep()


gert::git_add(".")
gert::git_commit("Initial commit")

gert::git_commit("CPB-BIO-DAD-preprocessing-v1")
gert::git_commit("CPB-BIO-DAD-preprocessing-v2")

usethis::use_git_config(user.name = "Tu Nombre", user.email = "tu@email.com")

gert::git_branch()

usethis::git_default_branch_configure()   # renombra master a main si hace falta
usethis::use_github()

