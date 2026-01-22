FROM julia:1.10

WORKDIR /app
ENV JULIA_PROJECT=/app
ENV JULIA_DEPOT_PATH=/app/.julia_depot

COPY Project.toml Manifest.toml ./
RUN mkdir -p /app/.julia_depot && chmod -R 777 /app/.julia_depot
RUN julia --project=. -e 'using Pkg; Pkg.instantiate(; allow_autoprecomp=false)'

COPY . .

RUN julia --project=. -e 'using Pkg; Pkg.resolve(); Pkg.instantiate(; allow_autoprecomp=false)'

ENV JULIA_NUM_THREADS=2
ENV JULIA_PKG_PRECOMPILE_AUTO=0
EXPOSE 8000

CMD ["bash", "/app/entrypoint.sh"]
