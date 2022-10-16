FROM alpine:latest

RUN apk add --no-cache bash

WORKDIR /app

COPY target/x86_64-unknown-linux-musl/release/sage ./sage

ENV PATH="/app:$PATH"
