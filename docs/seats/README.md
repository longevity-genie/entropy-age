# Seats

The standing roles agent sessions take in this repository, run by the seats kit
(agent-seat-blueprints; `docs/DECISIONS.md` D3, D5). The seat table is `.blueprints.toml` at the root; a
session takes a seat with `/seats:seat NAME`, one session per seat.

| Seat | Role (the kit's) | Owns | Priming | Notes |
|---|---|---|---|---|
| `lead` | `lead` | `docs/ROADMAP.md`, the code, figures and research docs | [LEAD.md](LEAD.md) | [LEAD_NOTES.md](LEAD_NOTES.md) |
| `runner` | `runner` | `docs/RUN_STATUS.md`, `docs/RUN_LOG.md`, its cluster jobs | [RUNNER.md](RUNNER.md) | [RUNNER_NOTES.md](RUNNER_NOTES.md) |
