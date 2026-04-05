---
name: download-large-file
description: Download large files robustly with resume and retries. Prefer aria2c first, then fall back to wget on failure or when aria2c is unavailable.
argument-hint: [url] [output-path]
---

Download large files with a strict fallback policy.

1. Validate inputs and destination.
- Require a full download URL and target output path.
- Ensure parent directory exists: `mkdir -p "$(dirname "<output-path>")"`.
- If a partial file exists, keep it for resume.

2. Try `aria2c` first.
- Check availability: `command -v aria2c`.
- Use resumable multi-connection download:
```bash
aria2c \
  --continue=true \
  --max-connection-per-server=16 \
  --split=16 \
  --min-split-size=1M \
  --retry-wait=5 \
  --max-tries=10 \
  --timeout=60 \
  --dir "$(dirname "<output-path>")" \
  --out "$(basename "<output-path>")" \
  "<url>"
```

3. Fall back to `wget` if `aria2c` is missing or exits non-zero.
- Preferred one-shot command (automatic fallback):
```bash
if command -v aria2c >/dev/null 2>&1 && aria2c \
  --continue=true \
  --max-connection-per-server=16 \
  --split=16 \
  --min-split-size=1M \
  --retry-wait=5 \
  --max-tries=10 \
  --timeout=60 \
  --dir "$(dirname "<output-path>")" \
  --out "$(basename "<output-path>")" \
  "<url>"; then
  echo "download succeeded via aria2c"
else
  wget \
    --continue \
    --tries=10 \
    --waitretry=5 \
    --read-timeout=60 \
    --timeout=60 \
    --output-document="<output-path>" \
    "<url>"
fi
```

4. Verify success before reporting done.
- Confirm output file exists and size is greater than zero.
- If a checksum is provided, verify it (`sha256sum` preferred).
- Report which tool actually succeeded (`aria2c` or `wget`).

5. Failure reporting format.
- Include URL, output path, tool attempted, exit code, and last stderr lines.
- Suggest one retry with lower `aria2c` concurrency (for example `--split=4`).
