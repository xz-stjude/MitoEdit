FROM continuumio/miniconda3:latest

# Set working directory
WORKDIR /app

# Copy environment files
COPY web-environment.yml .
COPY talen-environment.yml .

# Create environments using conda
RUN conda env create -f web-environment.yml && \
    conda env create -f talen-environment.yml && \
    conda clean --all --yes

# Initialize conda in bash
RUN conda init bash && \
    echo "conda activate mitoedit-web" >> ~/.bashrc

# Expose port 80
EXPOSE 80

# Set environment variables
ENV PORT=80
# MITOEDIT_PASSWORD is required but not set here for security reasons
# When running the container, set it with: -e MITOEDIT_PASSWORD=your_secure_password

# Copy application code
COPY . .

# Start the web server
CMD ["conda", "run", "--no-capture-output", "-n", "mitoedit-web", "python", "web/main.py"]
