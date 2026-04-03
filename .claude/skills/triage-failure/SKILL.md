---
name: triage-failure
description: Triage any failed run, test, or command quickly and produce a root-cause hypothesis plus a concrete recovery command. Use when logs show errors, flaky behavior, or non-zero exits.
argument-hint: [failed-command-or-log-path]
---

Triage failures with speed and evidence.

1. Reproduce or inspect the exact failure first.
- Re-run the failing command if safe and fast.
- Otherwise inspect the provided logs and stderr text.

2. Classify the failure.
- Environment/configuration (missing dependency, wrong path, permissions).
- Data/input contract (invalid format, null/empty data, missing columns/files).
- Logic/bug (exception from project code).
- Resource/runtime (OOM, timeout, deadlock, excessive parallelism).
- External dependency/network/API failure.

3. Isolate minimal cause.
- Identify first failing boundary (command, function, module, file).
- Extract the smallest reproducible command or test.
- Separate symptom from root-cause candidate.

4. Produce remediation.
- Provide one primary fix and one fallback fix.
- Give exact command(s) to recover.
- Include prevention action (test, guard, retry/backoff, input validation, docs update).

5. Report in strict format.
- Failure signature.
- Root-cause hypothesis (confidence high/medium/low).
- Evidence inspected.
- Recovery command.
- Prevention action.
