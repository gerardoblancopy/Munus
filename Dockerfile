FROM julia:1.10

WORKDIR /app

COPY Project.toml Manifest.toml ./
RUN julia -e 'using Pkg; Pkg.instantiate()'

COPY . .

ENV JULIA_NUM_THREADS=2
EXPOSE 8000

CMD ["julia", "--project=.", "server.jl"]
