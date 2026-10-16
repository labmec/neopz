---
type: finding
status: reviewed
updated: 2026-07-02
confidence: high
classification: confirmed-bug
severity: minor
evidence-commit: 6ffd38b12
tags:
  - divfreebubbles
  - application
  - benchmarking
---

# App repo: peak-memory reporting is platform-dependent (bytes vs KB), distorting benchmark tables 1024× on macOS

## Repository evidence
`divfreebubbles/targets/iter_elast.cpp:44-50` [repo]: `getPeakMemoryMB()` returns `usage.ru_maxrss / 1024.0` with comment "ru_maxrss is in KB on Linux". Same helper pattern in dupl_connects. Values land in the last column of `results_*_memory_time.txt` and in stdout as "Memory usage: … MB".

## Runtime evidence
[run @ 852a5116c(+), 2026-07-02]: on this Mac, idiv=50 prints "Memory usage: 296304 MB" (≈289 real MB — i.e. the printed number is KiB); idiv=400 prints 1.01205e7 (≈9.65 GB real). The user's own `results_Elastic2D_memory_time.txt` (written Jul 2) shows the same magnitudes → their local benchmark data carries the same unit.

## Reference evidence
POSIX leaves `ru_maxrss` units unspecified; Linux reports **kilobytes** (`getrusage(2)`), macOS/BSD reports **bytes** (`getrusage(2)` BSD man page). Well-known portability trap.

## Assessment
Classification: **confirmed bug — application repo** (benchmark instrumentation, not NeoPZ). Severity minor for library correctness, but *material for research outputs*: memory columns produced on macOS overstate 1024×; cross-machine comparisons (commit "adjusting target to get the peak memory"/"run in dell" suggests mixed Linux/Mac usage) are inconsistent. Fix: `#ifdef __APPLE__ /1024/1024 else /1024`.

## Related
[[flow-iter-elast]] · [[flow-dupl-connects]]
