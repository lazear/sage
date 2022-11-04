FROM debian:bullseye-slim

RUN apt-get update && \
    apt-get install -y --no-install-recommends procps \
    && rm -rf /var/lib/apt /var/lib/dpkg /var/lib/cache /var/lib/log

WORKDIR /app

COPY target/x86_64-unknown-linux-gnu/release/sage /app/sage

ENV PATH="/app:$PATH"
