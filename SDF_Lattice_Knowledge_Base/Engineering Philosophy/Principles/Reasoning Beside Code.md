---
title: Reasoning Beside Code
tags: [philosophy, principle, documentation]
---

# Reasoning Beside Code

> When producing code, also produce the **surrounding reasoning**.

Code without reasoning is a fact without an argument. The fact may be correct; the absence of the argument means no one can verify that, defend it under refactoring, or extend it with confidence.

## What "surrounding reasoning" means

For each non-trivial piece of work, deliver:

- **The architecture.** What the components are, how they fit together, why this decomposition. Not in exhaustive detail — enough that the next reader does not have to reverse-engineer the layout.
- **The invariants.** What properties the system maintains. Where each is established. Where each is relied on.
- **The likely failure modes.** What could go wrong. What would it look like. What would the system do.
- **The testing strategy.** Especially the **adversarial** part — see [[Tests Should Make Programs Fail]] and [[Output Format]]. Which failure classes the tests target, and which they intentionally don't.

This is the work of [[Output Format]]'s four-section template, applied at any scale: a one-line edit doesn't need a treatise; a new module deserves a paragraph or three; a new subsystem deserves the full structure.

## Why this matters

### 1. It externalizes the engineer's mental model

The hardest part of joining a project late is rebuilding the original engineer's mental model. Recorded reasoning short-circuits that. Without it, every new contributor pays the same archaeology tax.

### 2. It makes the design auditable

A claim like "this can never deadlock" is unauditable on its own. The same claim with a brief argument — "no two locks are ever acquired in different orders, because all multi-lock acquisitions go through `acquire_pair()`, which orders them by address" — can be checked. Future changes that violate the argument become visible as violations.

### 3. It catches design flaws in the writing

The act of writing the reasoning often reveals the holes. A correctness argument you cannot complete is a correctness argument that **does not exist yet**, even if the code happens to work. Better to find that out at design time than at incident time.

### 4. It protects against future "small changes"

Most catastrophic regressions come from "small changes" made without understanding the load-bearing assumptions. Recorded reasoning makes assumptions visible. The change either respects them, updates them deliberately, or visibly violates them.

## Where the reasoning lives

Different scales of reasoning belong in different places:

- **System-wide invariants** → vault notes (this knowledge base).
- **Subsystem / module-level reasoning** → module-level doc comments and design docs.
- **Function-level invariants** → `///` doc comments. See [[Documentation as Truth]].
- **Single-line subtleties** → a brief code comment, only when the *why* would surprise a reader.

Match the scale of the reasoning to the scale of the decision. Treatises about trivial helpers are noise; a one-line comment on a load-bearing invariant is gold.

## What "reasoning" is NOT

- It is not "explain what the code does." Well-named code does that. The reasoning is *why* the code is shaped this way.
- It is not "narrate every function call." Most calls speak for themselves.
- It is not retroactive marketing. Documentation written after the fact to make a design *sound* coherent without actually verifying it is worse than no documentation at all — it lies.

## A small example

### Without reasoning

```rust
pub fn submit(&mut self, job: Job) -> Result<JobId, SubmitError> {
    let id = self.next_id;
    self.next_id += 1;
    self.queue.push_back((id, job));
    Ok(JobId(id))
}
```

### With reasoning

```rust
/// Submit a job to the queue, returning a unique [`JobId`] for tracking.
///
/// # Invariant
/// `JobId`s are strictly monotonically increasing within a single
/// `Scheduler` instance. Callers may rely on `id_a < id_b` implying
/// `job_a` was submitted before `job_b`.
///
/// # Errors
/// Returns [`SubmitError::QueueFull`] when the queue has reached
/// `MAX_PENDING_JOBS`; the caller may retry after `poll`.
///
/// # Note
/// Wraparound of the underlying `u64` counter is not handled — at
/// the documented submission rate this is centuries away. If the
/// rate increases, revisit.
pub fn submit(&mut self, job: Job) -> Result<JobId, SubmitError> {
    if self.queue.len() >= MAX_PENDING_JOBS {
        return Err(SubmitError::QueueFull);
    }
    let id = self.next_id;
    self.next_id += 1;
    self.queue.push_back((id, job));
    Ok(JobId(id))
}
```

Same code, much more legible argument. Future maintainers know what they are allowed to rely on, what would break the invariant, and which assumption (`u64` wraparound) is parked deliberately.

## Related

- [[Human Understanding First]]
- [[Reason About Correctness]]
- [[Output Format]]
- [[Documentation as Truth]]
