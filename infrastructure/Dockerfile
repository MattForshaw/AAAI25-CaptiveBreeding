# A pinned environment for reproducing CoBreeder runs.
#
# Build and run from the repository root:
#   docker build -f infrastructure/Dockerfile -t cobreeder .
#   docker run --rm cobreeder run data/dir/individuals.csv \
#       data/dir/pr-scaled.csv data/dir/corrals-a.csv \
#       ALL_PAIRS_PR_MIN_SQUARED demo 0 0 10
FROM python:3.11-slim

# Dependencies come from poetry.lock, so the image and a local `poetry install` resolve
# to identical versions. Installing into the system interpreter keeps the image small;
# there is nothing else in it to conflict with.
RUN pip install --no-cache-dir "poetry==1.8.5" \
 && poetry config virtualenvs.create false

WORKDIR /app
COPY pyproject.toml poetry.lock ./

# --only main omits asv, which is a benchmarking tool and not needed to solve.
# --no-root installs the dependencies without the project, which is copied next.
RUN poetry install --only main --no-root

COPY cobreeder/ ./cobreeder/
COPY data/ ./data/
COPY examples/ ./examples/
COPY experiments/ ./experiments/
COPY README.md LICENSE ./

# Records what the image actually resolved, so a run can be traced to an environment
# even if the image is rebuilt later.
RUN python - <<'PY' > /app/environment.txt
import sys, ortools, pandas
print("python  ", sys.version.split()[0])
print("ortools ", ortools.__version__)
print("pandas  ", pandas.__version__)
PY

ENTRYPOINT ["python", "cobreeder/solver.py"]
CMD ["--help"]
