FROM rust
WORKDIR /app
COPY . /app
RUN cargo run --release tests/config.json
ENTRYPOINT ["/app/target/release/sage"]
