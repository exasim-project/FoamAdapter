#!/usr/bin/env bash
# --------------------------------------------------------------------------------------------
# Monitor the status of a LRZ GitLab CI pipeline on TUM COMA cluster for a given pipeline ID.
#
# Uses environment variables defined in the GitHub Actions job:
#   LRZ_GROUP, LRZ_HOST, REPO_NAME, MAX_WAIT_MINUTES, LRZ_GITLAB_PROJECT_TOKEN
#
# Usage:
#   .github/scripts/wait_pipeline.sh <pipeline_id>
# --------------------------------------------------------------------------------------------

set -euo pipefail

if [ $# -lt 1 ]; then
  echo "Usage: $0 <pipeline_id>"
  exit 1
fi

PIPELINE_ID=$1
GROUP="${LRZ_GROUP:?LRZ_GROUP not set}"
PROJECT="${REPO_NAME:?REPO_NAME not set}"
PROJECT_TOKEN="${LRZ_GITLAB_PROJECT_TOKEN:?LRZ_GITLAB_PROJECT_TOKEN not set}"
HOST="${LRZ_HOST:?LRZ_HOST not set}"
WAIT_MINUTES="${MAX_WAIT_MINUTES:-60}"

pipeline_url="https://${HOST}/${GROUP}/${PROJECT}/-/pipelines/${PIPELINE_ID}"
echo "Monitoring LRZ GitLab CI pipeline: $pipeline_url"

for i in $(seq 1 "$WAIT_MINUTES"); do
  status=$(curl -s \
    --header "PRIVATE-TOKEN: ${PROJECT_TOKEN}" \
    "https://${HOST}/api/v4/projects/${GROUP}%2F${PROJECT}/pipelines/${PIPELINE_ID}" \
    | jq -r '.status')

  echo "[$i] $PROJECT pipeline status: $status"

  case "$status" in
    success)
      echo "$PROJECT CI pipeline succeeded"
      exit 0
      ;;
    failed|canceled|skipped)
      echo "$PROJECT CI pipeline finished with status: $status"
      exit 1
      ;;
  esac

  sleep 60
done

echo "Timed out after $WAIT_MINUTES minutes waiting for $PROJECT CI pipeline"
exit 1
