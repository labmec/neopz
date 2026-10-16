---
type: finding
status: reviewed
updated: 2026-07-02
confidence: high
classification: confirmed-bug
severity: major
evidence-commit: 6ffd38b12
tags:
  - neopz
  - cpp
  - error-handling
---

# DebugStop() throws a messageless std::bad_exception unconditionally — including Release builds

## Repository evidence (verified first-hand)
`Common/pzerror.cpp:10-29`: `DebugStopImpl` prints "Your chance to put a breakpoint at file:line" to `PZError` and then `throw std::bad_exception();` — the `#ifdef PZDEBUG` guards are **commented out** (:15,:17), so this runs in every configuration. Call-site scale: ~3,029 `DebugStop()` occurrences vs ~51 `throw` and ~9 `catch` library-wide [agent counts].

## Why it matters
- The library's universal assertion mechanism behaves as an unrecoverable, message-free exception in production: uncaught (9 catches for 3k sites), it terminates with no file:line in the exception object (only on cerr, which GUIs/batch systems may swallow).
- `std::bad_exception` is semantically wrong (it exists for exception-specification violations), and carrying no payload makes programmatic handling impossible — downstream apps (e.g. divfreebubbles) cannot distinguish "hanging-node unsupported here" from "matrix singular".
- Interacts with [[finding-build-config-gaps]]: in the default RelWithDebInfo build, the `#ifdef PZDEBUG` *pre-checks* are compiled out while unguarded `DebugStop()`s still throw — an inconsistent middle state.

## Classification
**Confirmed defect (accidental complexity)** — not a domain requirement. The commented-out guard shows the release behavior is unintentional or at least unreviewed.

## Suggested improvement / risk
Introduce a typed `TPZFatalError{file,line,msg}` (or abort in debug, throw typed in release); mechanical sed-scale replacement, low risk; enables meaningful catches at analysis boundaries. Consider keeping `DebugStop` name as macro alias for compatibility.

## Related
[[finding-build-config-gaps]] · [[material-system]] · CPP_TECHNICAL_REVIEW §2-H1
