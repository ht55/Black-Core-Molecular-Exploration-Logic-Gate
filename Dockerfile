# --- Stage 1: Frontend Build ---
FROM node:18 AS build-stage
WORKDIR /app/frontend
COPY frontend/package*.json ./
RUN npm install --legacy-peer-deps
COPY frontend/ ./
RUN npm run build

# --- Stage 2: Backend & Runtime ---
FROM python:3.10-slim

RUN useradd -m -u 1000 user
USER user
ENV HOME=/home/user \
    PATH=/home/user/.local/bin:$PATH

WORKDIR $HOME/app

USER root
RUN apt-get update && apt-get install -y \
    libxrender1 libxext6 libsm6 \
    && rm -rf /var/lib/apt/lists/*
USER user

COPY --chown=user:user backend/requirements.txt .
RUN pip install --no-cache-dir --user -r requirements.txt

COPY --from=build-stage --chown=user:user /app/frontend/build ./build
COPY --chown=user:user backend/ .

EXPOSE 7860


CMD ["uvicorn", "main:app", "--host", "0.0.0.0", "--port", "7860"]