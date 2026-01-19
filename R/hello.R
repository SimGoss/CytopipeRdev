# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

hello <- function() {
  print("Hello, world!")
}

usethis::use_git_config(user.name = "Cyril-Etienne", user.email = "cyril.etienne@cea.fr")

## create a personal access token for authentication:
usethis::create_github_token()
## in case usethis version < 2.0.0: usethis::browse_github_token() (or even better: update usethis!)

## Note for Linux users:
## credentials::set_github_pat() (in line 34) might store your PAT in a memory cache that
## expires after 15 minutes or when the computer is rebooted. You thus may wish to do
## extend the cache timeout to match the PAT validity period:
usethis::use_git_config(credential.helper="cache --timeout=2600000") #> cache timeout ~30 days

## set personal access token:
credentials::set_github_pat("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")

credentials::set_github_pat("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")



