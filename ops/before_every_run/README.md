# Before Every Run

This folder stores the shared run-memory for remote execution on `ubuntu-tail`.

Both direct remote operation and Mac-led SSH orchestration are valid. This
remote folder remains the canonical run-memory surface either way; the Mac
mirror is for coordination and handoff.

## Rule

Before launching a new remote run, read:

- `LATEST.md`
- at least the newest relevant journal entry

Prefer reading the full journal when the task belongs to the same workflow family.

After the run finishes or fails:

- update the journal
- refresh `LATEST.md`

## Structure

- `LATEST.md`
- `journal/`
