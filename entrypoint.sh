#!/usr/bin/env bash
set -euo pipefail

export JULIA_PROJECT=/app
export JULIA_DEPOT_PATH=/app/.julia_depot
export JULIA_PKG_PRECOMPILE_AUTO=0

mkdir -p /app/.julia_depot

exec julia --compiled-modules=no --project=. /app/server.jl
