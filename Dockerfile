FROM alpine:latest

WORKDIR /app

COPY target/x86_64-unknown-linux-musl/release/sage ./sage

ENV PATH="/app:$PATH"
