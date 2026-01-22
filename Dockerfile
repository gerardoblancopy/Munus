FROM julia:1.10

WORKDIR /app
ENV JULIA_PROJECT=/app

COPY Project.toml Manifest.toml ./
RUN julia --project=. -e 'using Pkg; Pkg.instantiate(; allow_autoprecomp=false)'

COPY . .

ENV JULIA_NUM_THREADS=2
ENV JULIA_PKG_PRECOMPILE_AUTO=0
ENV JULIA_COMPILED_MODULES=0
EXPOSE 8000

CMD ["julia", "--compiled-modules=no", "--project=.", "server.jl"]
