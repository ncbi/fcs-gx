# used for passing git-info to main.cpp via compiler flag in src/CMakeLists.txt
# GP-32810: If branch is production, suppress git-revision (--abbrev=0), and no --dirty

branch_name=$(git rev-parse --abbrev-ref HEAD)

if [ "$branch_name" == "production" ]; then
    # e.g. v0.4.0
    str=$(git describe --always --tags --abbrev=0)
elif [ "$branch_name" == "develop" ]; then
    # e.g. v0.4.0-211-g2ee769ba-dirty
    # https://semver.org/#spec-item-9
    str=$(git describe --always --tags --dirty)
else
    # For feature-branches, capture the branch-name:
    # e.g. v0.4.0-211-g2ee769ba-dirty+branch--task-capture-branch-name-in-git-describe-sh
    # https://semver.org/#spec-item-10
    alnum_branch_name=$(echo "$branch_name" | sed 's/[^a-zA-Z0-9]/-/g')
    str=$(git describe --always --tags --dirty)"+branch--${alnum_branch_name}"
fi

printf "\"$str\"" # output must be quoted
