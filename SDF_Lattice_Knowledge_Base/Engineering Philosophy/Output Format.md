---
title: Output Format
tags: [philosophy, output, template]
---

# Output Format

For non-trivial design or implementation work, structure the deliverable in **four sections**, in this order. The format exists because each section catches a class of failure the others cannot.

## 1. Correctness model

State, in plain language, **what the system must do for it to be considered correct**. Not what it must compute — what property must hold. Examples:

- "For any valid lattice description, the solver produces a vertex set whose neighborhood graph is connected."
- "For any GPU output, there exists a CPU computation in `core` that produces the same answer within a documented tolerance."
- "Calling `flush()` on an open writer makes all prior `write()` calls observable to subsequent readers, even after process death."

The correctness model is the single contract the rest of the work is justified against. If you can't write it in one paragraph, the design is not yet understood.

## 2. Invariants and assumptions

List the **load-bearing claims** the design relies on. Two flavors:

- **Invariants the system maintains** — properties true at every observation point (e.g., "the index map is always a bijection", "no `Port` value is ever zero").
- **Assumptions about the environment** — properties relied on but not enforced by the code (e.g., "the GPU adapter supports compute pipelines", "the input file is no larger than 4 GiB").

Each invariant should name **where it is established** (constructor, validation step, type definition) and **where it is relied on** (callers, downstream functions).

This section is doing the work of [[Reason About Correctness]] — externalizing the argument so a reader can audit it.

## 3. Design and implementation

The actual code, plus the prose that lets a reader follow it. Conventions:

- **Why before how.** Each major component opens with one to three sentences of intent.
- **Narrate the data movement.** A reader should be able to trace a value from input to output without running the code.
- **Surface the type-system arguments.** When a type makes an invalid state impossible, say so. The compiler enforcing it is not the same as a reader knowing it.
- **Keep cleverness on a leash.** Cleverness is fine; undocumented cleverness is debt.

This section embodies [[Human Understanding First]] and [[Explicit Not Magical]].

## 4. Test strategy

Split into four labeled subsections. Each test class targets a different failure mode; omitting any of them leaves a class uncovered.

### a. Ordinary behavioral tests
- Confirm the design works on representative valid inputs.
- Catch outright wrong implementations.
- These are necessary but not sufficient.

### b. Edge and boundary tests
- Empty, full, single-element, max, max+1, zero, negative, exact equality at thresholds.
- Catch off-by-one, fence-post errors, mishandled empty cases.
- See [[Edge Cases and Properties]].

### c. Adversarial / torture tests
- Combine edges unpleasantly (empty + max, malformed + valid mixed).
- Repeated operations, idempotence violations, unusual orderings.
- Contradictory options. Invalid state transitions attempted from outside.
- Tests that **would embarrass a weak design**.
- Catch the bugs ordinary usage will never reveal — the entire point of [[Normal Usage Is Not Testing]].

### d. Regression tests
- Each historical bug becomes a permanent test.
- Even if the bug seems unlikely to recur — *especially* if it does.
- Tests preserve hard-won knowledge that is otherwise only in commit messages and people's heads.
- See [[Regression Discipline]].

For **each class**, briefly state **what failure mode it is intended to expose**. A test list without that annotation is a list of activities, not a defense.

## When to skip the format

For trivial requests — a one-line edit, a question, a renamed variable — don't ceremoniously apply the four sections. But still:

- Favor [[Sharp Oracles|sharp oracles]] in any tests written.
- Surface invariants when they're load-bearing.
- Don't smuggle in cleverness without a comment.

The format scales with the stakes. The discipline doesn't.

## Related

- [[Engineering Philosophy MOC]]
- [[Reasoning Beside Code]]
- [[Sharp Oracles]]
- [[Three Levels of Tests]]
