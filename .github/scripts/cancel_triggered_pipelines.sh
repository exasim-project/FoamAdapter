#!/usr/bin/env bash
# ---------------------------------------------------------------------------------------
# Cancel running/pending LRZ GitLab CI pipelines on TUM COMA cluster for a given branch.
# Only cancel if the pipeline does not have the NEON_BRANCH variable set.
#
# Uses environment variables defined in the GitHub Actions job:
#   LRZ_GROUP, LRZ_HOST, REPO_NAME, LRZ_GITLAB_PROJECT_TOKEN
#
# Usage:
#   .github/scripts/cancel_pipelines.sh <branch>
# ---------------------------------------------------------------------------------------

set -euo pipefail

if [ $# -lt 1 ]; then
  echo "Usage: $0 <branch>"
  exit 1
fi

BRANCH=$1
GROUP="${LRZ_GROUP:?LRZ_GROUP not set}"
PROJECT="${REPO_NAME:?REPO_NAME not set}"
TOKEN="${LRZ_GITLAB_PROJECT_TOKEN:?LRZ_GITLAB_PROJECT_TOKEN not set}"
HOST="${LRZ_HOST:?LRZ_HOST not set}"

echo "Checking running/pending CI pipelines for branch: $BRANCH in project: $GROUP/$PROJECT"

# Fetch pipelines for the branch
response=$(curl -s -w "%{http_code}" -o response.json \
  --header "PRIVATE-TOKEN: ${TOKEN}" \
  "https://${HOST}/api/v4/projects/${GROUP}%2F${PROJECT}/pipelines?ref=${BRANCH}&order_by=id&sort=desc")

http_code="${response:(-3)}"

if [[ "$http_code" != "200" ]]; then
  echo "GitLab API request failed with HTTP status $http_code"
  cat response.json
  exit 1
fi

# Ensure response is a JSON array
if ! jq -e 'type=="array"' response.json >/dev/null 2>&1; then
  echo "Unexpected response from GitLab API (not a JSON array)"
  cat response.json
  exit 1
fi

# Extract running/pending pipeline IDs
pipeline_ids=$(jq -r '.[] | select(.status=="running" or .status=="pending") | .id' response.json)
if [ -z "$pipeline_ids" ]; then
  echo "No running/pending CI pipelines to check"
else
  for id in $pipeline_ids; do
    echo "Checking pipeline $id for NEON_BRANCH variable"
    vars=$(curl -s \
      --header "PRIVATE-TOKEN: ${TOKEN}" \
      "https://${HOST}/api/v4/projects/${GROUP}%2F${PROJECT}/pipelines/$id/variables")

    neon_branch=$(echo "$vars" | jq -r '.[] | select(.key=="NEON_BRANCH") | .value' || true)

    if [ -z "$neon_branch" ]; then
      echo "Canceling pipeline $id (NEON_BRANCH is null)"
      curl -s --request POST \
        --header "PRIVATE-TOKEN: ${TOKEN}" \
        "https://${HOST}/api/v4/projects/${GROUP}%2F${PROJECT}/pipelines/$id/cancel" >/dev/null
    else
      echo "Keeping pipeline $id (NEON_BRANCH=$neon_branch)"
    fi
  done
fi

echo "Done."
