# Decisions

The calls settled with the maintainer, one section each, `## Dn — the call (date)`, newest last: what
was decided, why, and the alternative it rejected. A call is never edited after the fact; a later one
that changes it says which it amends.

## D1 — Push `main` directly, no pull requests (2026-11-09)

The maintainer works alone here: sessions commit and push `main` as a fast-forward without asking,
and open a pull request only when asked or when the push is refused. The seats kit's adoption
wizard lands its files as a pull request whose merge is the approval; that was rejected for the
adoption commit too, in favour of the maintainer's standing rule. Force pushes and rewriting pushed
history stay forbidden.

## D2 — The existing lowercase docs keep their names (2026-11-09)

`CLAUDE.md` names new files under `docs/` in UPPER_SNAKE. `docs/paper.md`, `docs/model.md` and
`docs/extension.md` predate that rule and the README links them, so they keep their names; renaming
them was rejected as churn that breaks links. A default taken at adoption, listed for review.

## D3 — Seats: one lead, tracked, on a new Zulip organisation (2026-11-09)

The repository adopts the seats kit (agent-seat-blueprints, v0.1.144) with one seat, `lead`: it owns
`docs/ROADMAP.md` and builds it with the maintainer. Rejected for now: a liaison following the
upstream repository and GEO, a reviser for a manuscript, a presenter; each can be added later.

- **Tracked, not an overlay:** `.claude/settings.json` pins the kit and `.blueprints.toml` is
  committed, so every clone gets the plugin. The overlay (nothing of the kit tracked) was rejected,
  although the kit's repository is private and this one is public: the address is named here on the
  maintainer's call, and a clone without access to the kit simply does not load the plugin.
- **Messenger:** a new Zulip organisation of the maintainer's own, with its secretary seat and the
  administrator key on another machine than this repository's. Rejected: no messenger (seats on one
  machine only), and someone else's organisation.
- **Desk:** a new desk, `longevity-genie/desk` (private), the organisation's, made on the
  secretary's machine. Rejected: no desk (questions waiting in the chat).
- **The lead's grant:** commit, push `main` as a fast-forward, and release tags. Rejected: commit
  only (the kit's safest default), and push without tags.
- **The leak gate** guards this machine's hostname and home path, and the Zulip server's address
  once it exists; not the maintainer's login, which `pyproject.toml` already carries.

## D4 — The research programme: reimplement, add flavours, validate on contrasts (2026-11-09)

The maintainer (typos corrected): "we are inspired by major ideas of the stochastic clock ... So we start from
implementing what is there in a better way and then we capitalize on the idea to add new flavour of
stochastic clock in transcriptomics and methylation where we will try to explore stochasticity on
selected aspects (hallmark genes) and selected cell types (for example exercise can select some cells
that have better fitness on some genes) and then we will try to validate them to see if there is a
signal in samples with quite contrast (for example exercised vs lazy of known ages)."

The sources are Meyer & Schumacher (2024) and David H. Meyer's EHA 2026 talk (BiT age, stochastic
clocks from one ground-state sample, interventions, reprogramming, diapause, NeuronAge). The order is
binding: a new flavour builds on a reproduction that has been checked against the original. Rejected:
a pure reproduction with no extension, and starting from the extension before the reproduction
stands. The lead's roadmap is drafted in these three steps (`docs/seats/LEAD.md`).
