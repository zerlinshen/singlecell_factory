---
name: orchestrate-large-task
description: Orchestrate large, cross-domain tasks by splitting work into phases and delegating independent work packages to typed subagents in parallel when policy and user intent permit delegation. Use when a task spans multiple modules/files, mixes build+debug+test+docs, or has independent subtasks where parallel execution materially reduces turnaround time.
---

Run a large task as a controlled composition workflow.

1. Trigger only when orchestration is justified.
- Use when work is cross-domain or likely touches 3+ files / 2+ modules.
- Do not trigger for small linear tasks; route directly to one specialized skill.

2. Execute five explicit phases.
- `Intake`: restate objective, constraints, and definition of done.
- `Decompose`: split into independent work packages with owners, inputs, outputs.
- `Delegate`: assign packages to typed subagents (`Scout`, `Builder`, `Tester`, `Reviewer`).
- `Integrate`: merge outputs, resolve conflicts, normalize contract compatibility.
- `Verify`: run targeted checks and summarize residual risks.

3. Enforce boundaries and anti-conflict rules.
- Delegate only when runtime policy allows subagent usage.
- Never assign overlapping write scopes to concurrent subagents.
- Keep blocking critical-path tasks in the main agent.
- Timebox delegated work and require concise deliverables.
- Reclaim stalled tasks immediately in the main agent.

4. Apply composition routing to prevent over/under-use.
- Pipeline execution/failure handling: call `/execute-and-recover-pipeline`.
- Module development/integration/refactor/dependency updates: call `/develop-and-integrate-module`.
- Review/output validation/docs alignment: call `/review-and-validate-quality`.
- Performance optimization/regression gates: call `/optimize-and-guard-performance`.
- Large file transfer prerequisites: call `/download-large-file`.

5. Keep cadence and reporting standardized.
- Track each phase with `pending | in_progress | done | blocked`.
- For each subagent, record scope, changed files, verification results, open risks.
- Final report must include completed work, unresolved items, residual risks, and exact verification commands.
