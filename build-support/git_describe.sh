# used for passing git-info to main.cpp via compiler flag in src/CMakeLists.txt
# GP-32810: If branch is production, suppress git-revision (--abbrev=0), and no --dirty
str=$(git describe --always --tags $([ $(git rev-parse --abbrev-ref HEAD) = "production" ] && echo "--abbrev=0" || echo "--dirty"))
printf "\"$str\"" # output must be quoted
