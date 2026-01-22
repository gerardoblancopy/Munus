#!/usr/bin/env bash
set -euo pipefail

export JULIA_PROJECT=/app
export JULIA_DEPOT_PATH=/tmp/julia
export JULIA_PKG_PRECOMPILE_AUTO=0

mkdir -p /tmp/julia

exec julia --compiled-modules=no --project=. /app/server.jl
