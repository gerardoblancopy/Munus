FROM julia:1.10

WORKDIR /app
ENV JULIA_PROJECT=/app
ENV JULIA_DEPOT_PATH=/app/.julia_depot

# Ensure depot directory exists and is writable
RUN mkdir -p /app/.julia_depot && chmod -R 777 /app/.julia_depot

COPY . .

# Resolve, Instantiate, and Precompile in one go to ensure consistency
RUN julia --project=. -e 'using Pkg; Pkg.resolve(); Pkg.instantiate(; allow_autoprecomp=false); Pkg.precompile()'

ENV JULIA_NUM_THREADS=2
EXPOSE 8000

CMD ["julia", "--project=.", "server.jl"]
