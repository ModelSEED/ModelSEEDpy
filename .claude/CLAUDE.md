# Claude Code Universal Guidelines

This file sets baseline expectations for Claude Code across all Hall-lab repos.
Individual repos may append their own conventions in a section after this one.

## Core behavior

- **Be specific.** Include file paths with line numbers, function names, error messages. State exactly what you changed and why.
- **Be honest.** If something is incomplete, uncertain, or failed, say so. Don't paper over problems.
- **Don't over-engineer.** Write the code the task needs — no speculative abstractions, no unrequested refactors, no "while we're here" cleanups.
- **Read before writing.** Don't propose changes to code you haven't read. Understand the existing pattern before adding to or modifying it.
- **Match existing conventions.** Follow the style, naming, and architecture already present in the repo. If you deviate, explain why.
- **Confirm destructive actions.** For anything hard to reverse (file deletion, force push, database truncate, rewriting history, destructive migrations), stop and ask before proceeding.

## File organization

When a task produces artifacts (design docs, PRDs, audits, plans, research notes), save them under `agent-io/`:

```
agent-io/
├── prds/<prd-name>/    # Product requirements: humanprompt.md, fullprompt.md, data.json
├── audits/             # Audit reports (YYYY-MM-DD-named)
├── docs/               # Architecture and usage documentation
├── plans/              # Design and implementation plans
└── research/           # Research notes and summaries
```

Code files go in the repo's existing source tree. Do not create code under `agent-io/`.

## Git hygiene

- Commits represent one logical unit of work. Don't bundle unrelated changes.
- Never force-push to `main`/`master` without explicit user approval.
- Never skip hooks (`--no-verify`) or bypass signing unless explicitly asked.
- Stage files by name when possible; avoid broad `git add .` to keep secrets out.
- Create new commits rather than amending published ones.

## Working with repo state

- If this repo has a `state/` directory, treat its contents as authoritative and mutate it only through the repo's documented API (Python modules, CLIs). Do not hand-edit JSON/YAML/SQLite state files.
- State files may include JSON logs, SQLite databases, YAML registries, activity sidecars. The repo's README or `.claude/CLAUDE.md` will describe its conventions.
