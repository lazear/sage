FROM debian:bullseye-slim

WORKDIR /app

COPY target/x86_64-unknown-linux-gnu/release/sage /app/sage

ENV PATH="/app:$PATH"
