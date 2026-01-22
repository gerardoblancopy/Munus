FROM julia:1.10

WORKDIR /app
ENV JULIA_PROJECT=/app
ENV JULIA_DEPOT_PATH=/tmp/julia

COPY Project.toml Manifest.toml ./
RUN mkdir -p /tmp/julia && chmod -R 777 /tmp/julia
RUN julia --project=. -e 'using Pkg; Pkg.instantiate(; allow_autoprecomp=false)'

COPY . .

ENV JULIA_NUM_THREADS=2
ENV JULIA_PKG_PRECOMPILE_AUTO=0
EXPOSE 8000

CMD ["julia", "--project=.", "server.jl"]
