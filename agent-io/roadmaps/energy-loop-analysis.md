---
project: modelseedpy
arc: energy-loop-analysis
title: Energy Loop Detection, Resolution & Database Cleaning
total_estimate: 6w
created: 2026-06-19
updated: 2026-06-19
status: not_started
milestones:
  - id: m1a
    title: Get loop detection fully working
    estimate: 1w
    status: not_started
    parallel_group: g1
    started: null
    completed: null
    handoff: null
    progress: null
  - id: m1b
    title: Evaluate & apply AI directionality constraints
    estimate: 1w
    status: not_started
    parallel_group: g1
    started: null
    completed: null
    handoff: null
    progress: null
  - id: m2
    title: Characterize surviving loops
    estimate: 1w
    status: not_started
    parallel_group: null
    started: null
    completed: null
    handoff: null
    progress: null
  - id: m3
    title: Generate impact data for remaining loops
    estimate: 1w
    status: not_started
    parallel_group: null
    started: null
    completed: null
    handoff: null
    progress: null
  - id: m4
    title: Reasoning-model voting council
    estimate: 1w
    status: not_started
    parallel_group: null
    started: null
    completed: null
    handoff: null
    progress: null
  - id: m5
    title: Cumulative integration, re-evaluation, finalize
    estimate: 1w
    status: not_started
    parallel_group: null
    started: null
    completed: null
    handoff: null
    progress: null
---

# Energy Loop Detection, Resolution & Database Cleaning — Roadmap

**Context.** Models contain energy loops (thermodynamically impossible cycles). Today we
rediscover and eliminate them every build — unpredictable and order-dependent. We want to
remove the *right* reaction, and clean loops out of the database once so future models
don't inherit them. Two pieces exist: (a) loop-detection code (written, not tested);
(b) impact-evaluation code (assesses follow-on effects of changing a reaction: new dead
reactions, what else breaks, which phenotypes/compounds can no longer be
produced/degraded).

## Milestone 1A — Get loop detection fully working (parallel, 1w)

Harden/validate the detection code: discover multiple loops; for each, chase all variants
by knocking out every component and checking for alternative reactions that fill the slot;
eliminate each loop by adding constraints so it isn't rediscovered; checkpoint progress so
the long-running job can stop/restart without re-detecting known loops.

*Runs in parallel with M1B.*

## Milestone 1B — Evaluate & apply AI directionality constraints (parallel, 1w)

Thousands of AI-proposed directionality changes already exist (e.g. "you call this
reversible; I think it's irreversible"), not yet applied. Evaluate each change's impact
one by one, stack the least-impactful changes (re-evaluating as each is added), merge
changes that don't produce a large combined detrimental effect in the proper order while
tracking database impact and reviewing overly harmful ones. As loops are found in M1A,
cross-reference which loops these changes already resolve.

*Runs in parallel with M1A.*

## Milestone 2 — Characterize surviving loops (1w)

After the directionality pass, determine which loops remain and characterize what's left,
as input to a resolution strategy.

## Milestone 3 — Generate impact data for remaining loops (1w)

Run impact-evaluation on every reaction in each remaining loop, giving the voting council
its data (e.g. "chemically correct to remove, but removing it kills half the database").

## Milestone 4 — Reasoning-model voting council (1w)

Submit each remaining loop (thousands to tens of thousands; token budget adequate) to a
council of reasoning models with the M3 impact data; the council votes on the reaction
most likely incorrect — the one to remove to disrupt the loop. Run systematically.

## Milestone 5 — Cumulative integration, re-evaluation, finalize (1w)

Evaluate the cumulative effect of stacking the voted fixes (same method as the
directionality work); integrate fixes that hold up; for loops whose first choice was
rejected, re-run evaluation for alternatives. Result: a database with as few loops as
possible, ready to document and publish.
