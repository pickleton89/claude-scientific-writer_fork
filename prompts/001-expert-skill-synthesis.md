<objective>
Extract the expert structure of a skill domain: principles, judgments, decision rules, and failure modes that separate competent from excellent execution.

You are an expert practitioner and critical synthesizer. Your task is NOT to teach tutorials or list tools. Your task is to distill what experts internalize but rarely articulate.
</objective>

<input>
Determine the skill domain from one of these sources (in order of priority):

1. **User-provided domain**: If the user specifies a domain (e.g., "peer review", "data visualization"), use that
2. **SKILL.md context**: If a skill path is provided, read the SKILL.md to extract the domain name and scope

Example invocations:
- `/run-prompt 001 peer review`
- `/run-prompt 001 skills/plotting-libraries`
</input>

<research>
Use WebSearch to gather current best practices for this domain:
- Search for "[domain] best practices 2024" or "[domain] expert techniques"
- Search for "[domain] common mistakes" or "[domain] pitfalls"
- Search for "[domain] quality criteria"

Synthesize findings with your knowledge. Prioritize specificity over comprehensiveness.
</research>

<internal_analysis instruction="do not output this section">
Before writing output, analyze:

1. What are the TRUE goals of this skill? (not surface-level, the deep purpose)
2. What mental models do experts use?
3. What practices reliably improve outcomes?
4. What do novices and AI systems consistently get wrong?
5. How do experts judge quality in this domain?
</internal_analysis>

<output_structure>
Generate a single-page reference document (500 words max) with these sections:

## 1. Core Principles (3-5)
Non-negotiable fundamentals experts internalize. Avoid generic advice. Be specific to this domain.

## 2. Decision Framework
Key questions and trade-offs that guide expert decisions.
Include context-dependent rules using format: When [situation], do [action].
Use a table format for clarity.

## 3. Concrete Patterns (2-3)
Reusable, actionable patterns in this format:

**Pattern Name**
- [Situation] → [Action] → [Rationale]

## 4. Failure Modes
Common mistakes that degrade outcomes:
- Reasonable approaches that fail
- Over-engineering patterns
- Cargo-cult practices

Use a table with columns: Mistake | Why It Fails

## 5. Anti-Patterns
What should NEVER be done, even if it seems efficient or standard. Be specific.

## 6. Quality Signals
How experts recognize strong work:
- "Good enough" vs "excellent" distinction
- Subtle cues distinguishing amateur from professional output

Use a comparison table.

## Final Synthesis (Required)
Conclude with exactly this format:

**Do This:** [3-5 expert defaults as bullet points]

**Not That:** [3-5 high-impact mistakes to avoid]

**When In Doubt:** [One guiding heuristic - a single sentence]

**Quality Check:**
1. [Self-review question]
2. [Self-review question]
3. [Self-review question]
4. [Self-review question]
5. [Self-review question]
6. [Self-review question]
</output_structure>

<constraints>
- Prioritize specificity over completeness
- Prefer principles over tools
- State trade-offs explicitly
- Assume an intelligent, experienced reader
- No fluff or motivational language
- No generic advice that applies to everything
- Every point should be specific to this domain
</constraints>

<file_output>
Determine the skill folder path and save to:
`./skills/[skill-name]/references/expert_guide.md`

If the skill folder doesn't exist or can't be determined, save to:
`./analyses/[domain-name]-expert-synthesis.md`

Add this header to the file:
```markdown
# [Domain Name] - Expert Guide

> Principles, decision rules, and failure modes that separate competent from excellent execution.
> Generated via expert synthesis prompt.

---
```
</file_output>

<verification>
Before completing, verify:
- [ ] All 6 sections are present
- [ ] Final Synthesis includes all 4 required components
- [ ] Word count is under 500 (excluding headers)
- [ ] No generic advice - every point is domain-specific
- [ ] Concrete patterns use the Situation → Action → Rationale format
- [ ] File saved to correct location
</verification>
