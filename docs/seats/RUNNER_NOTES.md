# Runner notes

The state the runner's last holder left. State, not rules: the rules are the priming, `RUNNER.md` beside
this file, and the kit's role `classes/runner/RUNNER.md`.

## 2026-11-10, first seating

- Held by the same session as the lead, on antonkulaga's word.
- The cluster machine's `slurm.toml` written from the account's own rules; the share is in `RUNNER.md`.
- No job submitted, no job watcher armed. R1 (GSE179847, RM5) waits on the lead's layout for that series.

## 2026-11-11

- Every job goes in with `--comment 'seat=entropy-age/runner; what=<one line>; eta=<ISO-8601 or ?>; cuda=-'`, and
  `scontrol update JobId=N Comment=...` when the estimate moves (the account's rigger, for every seat on the account).
  `what` names no person and no private data: every login reads it.
- R1 done (GSE179847). R2 (the per-person score, RM6) running as job 6428234, comment set.
