#!/usr/bin/env bash

set -euo pipefail

usage() {
    cat <<'EOF'
Usage: assets/sh/xorf-submodule.sh [--remote]

Without flags, initialize or restore modules/xorf at the commit recorded by this
repository and keep only the top-level /src/ tree checked out.

With --remote, update modules/xorf to the latest commit from the branch declared
in .gitmodules and keep the same sparse checkout. The superproject will then see
modules/xorf as modified until you commit the new submodule pointer.
EOF
}

if [[ $# -gt 1 ]]; then
    usage >&2
    exit 1
fi

remote_mode=0
if [[ $# -eq 1 ]]; then
    case "$1" in
        --remote)
            remote_mode=1
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            usage >&2
            exit 1
            ;;
    esac
fi

repo_root="$(git rev-parse --show-toplevel)"
cd "$repo_root"

submodule_path="modules/xorf"
submodule_key="submodule.${submodule_path}"
branch="$(git config -f .gitmodules --get "${submodule_key}.branch" || printf 'master')"

git submodule sync -- "$submodule_path"

update_args=(submodule update --init)
if [[ $remote_mode -eq 1 ]]; then
    update_args+=(--remote --depth 1)
fi
update_args+=(-- "$submodule_path")

git "${update_args[@]}"
git -C "$submodule_path" sparse-checkout set --no-cone '/src/'

printf 'xorf commit: %s\n' "$(git -C "$submodule_path" rev-parse HEAD)"
printf 'xorf branch: %s\n' "$branch"
printf 'checked out paths:\n'
find "$submodule_path" -mindepth 1 -maxdepth 1 ! -name .git -printf '  %f\n' | sort
