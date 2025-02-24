FROM condaforge/mambaforge:latest

# Set working directory
WORKDIR /app

# Copy environment files
COPY web-environment.yml .
COPY talen-environment.yml .

# Create environments using mamba
RUN mamba env create -f web-environment.yml && \
    mamba env create -f talen-environment.yml && \
    mamba clean --all --yes

# Initialize conda in bash
RUN conda init bash && \
    echo "conda activate mitoedit-web" >> ~/.bashrc

# Expose port 80
EXPOSE 80

# Set environment variables
ENV PORT=80

# Copy application code
COPY . .

# Start the web server
CMD ["conda", "run", "--no-capture-output", "-n", "mitoedit-web", "python", "web/main.py"]
