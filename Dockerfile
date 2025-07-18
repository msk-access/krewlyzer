# syntax=docker/dockerfile:1
FROM python:3.10-slim

# Install uv (fast Python package/dependency manager)
RUN pip install --no-cache-dir uv

# Set workdir
WORKDIR /app

# Copy project files
COPY . /app

# Install dependencies and package in a uv venv
RUN uv venv .venv && \
    . .venv/bin/activate && \
    uv pip install .

# Set environment variables for uv venv activation
ENV PATH="/app/.venv/bin:$PATH"

# Default command: show CLI help
ENTRYPOINT ["krewlyzer"]
CMD ["--help"]
