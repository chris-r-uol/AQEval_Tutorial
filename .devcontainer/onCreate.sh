#!/usr/bin/env bash
# Finishes the workspace: the Python packages, and RStudio's starting folder.
#
# This is onCreateCommand, so a prebuild includes it. The R packages are not
# installed here: the r-packages feature in devcontainer.json does that while
# the image is built.
#
# Nothing here stops the build. A student with a working R and a broken Python,
# or the reverse, can still do the tutorial in the other language, and the
# summary at the end says plainly which is usable.
set -uo pipefail

workspace="$(pwd)"
venv="$HOME/.venvs/aqeval"

echo "==> Python packages"
python3 -m venv "$venv"
"$venv/bin/pip" install --no-cache-dir --quiet --upgrade pip
"$venv/bin/pip" install --no-cache-dir --quiet -r python/requirements.txt
# The notebook asks for a kernel named "python3". Registering this environment
# under that name means it is picked without the student having to choose.
"$venv/bin/python" -m ipykernel install --user --name python3 \
	--display-name "Python 3 (AQEval tutorial)" >/dev/null

echo "==> RStudio"
# RStudio Server opens in the home folder by default, which here is empty.
# Opening in the repository instead lets .Rprofile put the tutorial script in
# front of the student. The lines are replaced, not appended, so that running
# this again leaves one copy.
config=/etc/rstudio/rsession.conf
sudo touch "$config"
sudo sed -i '/^session-default-working-dir=/d;/^session-default-new-project-dir=/d' "$config"
printf 'session-default-working-dir=%s\nsession-default-new-project-dir=%s\n' "$workspace" "$workspace" |
	sudo tee -a "$config" >/dev/null
echo "    opens in $workspace"

echo "==> Environment"
"$venv/bin/python" - <<'PY'
import importlib
for name in ["aqeval", "pandas", "numpy", "pygam", "plotly", "matplotlib"]:
    try:
        print(f"    python  {name:12s} {importlib.import_module(name).__version__}")
    except Exception as error:
        print(f"    python  {name:12s} UNUSABLE ({error})")
PY
# Load each package rather than reading its version number: a package can be
# installed and still fail to load, as AQEval does without Java.
Rscript -e 'cat(sprintf("    R       %-12s %s\n", "base", getRversion()))
  for (p in c("openair", "AQEval", "dplyr", "rstudioapi", "worldmet")) {
    ok <- tryCatch({ suppressMessages(library(p, character.only = TRUE)); TRUE },
                   error = function(e) FALSE)
    v <- if (ok) as.character(packageVersion(p)) else "UNUSABLE"
    cat(sprintf("    R       %-12s %s\n", p, v))
  }' || true
