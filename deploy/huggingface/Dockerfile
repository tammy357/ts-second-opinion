# Hugging Face Space (Docker SDK) for the marimo demo.
FROM python:3.12-slim
RUN apt-get update && apt-get install -y --no-install-recommends git libxrender1 libxext6 && rm -rf /var/lib/apt/lists/*
RUN useradd -m -u 1000 user
USER user
ENV HOME=/home/user PATH=/home/user/.local/bin:$PATH
WORKDIR /home/user
RUN git clone --depth 1 https://github.com/tammy357/ts-second-opinion app && pip install --no-cache-dir --user -e ./app
WORKDIR /home/user/app
EXPOSE 7860
CMD ["marimo", "run", "apps/demo.py", "--host", "0.0.0.0", "--port", "7860", "--headless"]
