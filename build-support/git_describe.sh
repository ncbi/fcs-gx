# used for passing git-info to main.cpp via compiler flag in src/CMakeLists.txt
# GP-32810: If branch is production, suppress git-revision (--abbrev=0), and no --dirty

script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
branch_name=$(git rev-parse --abbrev-ref HEAD 2>/dev/null)

if [ ! -n "$branch_name" ] && [[ -f $script_dir/../VERSION ]]; then
    # Not a git-repo (git rev-parse failed), but have VERSION file.
    str=$(head -n1 $script_dir/../VERSION | tr -d '\n') 

elif [ ! -n "$branch_name" ]; then
    str="unknown"

elif [ "$branch_name" == "production" ] || [ "$branch_name" == "release" ]; then
    # NB: 'release' is the default branch in fcs-gx repo.
    #
    # e.g. v0.4.0
    str=$(git describe --always --tags --abbrev=0)

elif [ "$branch_name" == "develop" ] || [ "$branch_name" == "main" ] || [ "$branch_name" == "HEAD" ]; then
    # NB: 'HEAD' is the branch-name that ends up in the singularity build.
    #
    # e.g. v0.4.0-211-g2ee769ba-dirty
    # https://semver.org/#spec-item-9
    str=$(git describe --long --dirty)

else
    # For feature-branches, capture the branch-name:
    # e.g. v0.4.0-211-g2ee769ba-dirty+branch--task-capture-branch-name-in-git-describe-sh
    # https://semver.org/#spec-item-10
    alnum_branch_name=$(echo "$branch_name" | sed 's/[^a-zA-Z0-9]/-/g')
    str=$(git describe --long --dirty)"+branch--${alnum_branch_name}"
fi

printf "\"$str\"" # output must be quoted
