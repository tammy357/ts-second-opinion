---
description: Turn a plain-English chemistry rule into a checked, gated validator for TS guesses
argument-hint: "\"rule in plain English\""
---
You are adding a rule to TS Second Opinion using the ts-second-opinion MCP tools.

Rule: $ARGUMENTS

1. Call list_rules and get_rule_source on one existing rule (for example "chirality") to see the expected style.
2. Write check(rxn, xyz) -> (passed, detail) in Python for the rule. numpy is `np`, math is `math`; no imports.
   Name atoms as symbol + (index + 1), e.g. "C3", in the detail string.
3. Call gate_check. If it is not accepted, read the feedback and the examples of real TSs it rejected, revise the code,
   and call gate_check again. Stop after 4 attempts.
4. If accepted, call save_rule with a short snake_case name (use_in_triage=False unless the user asked otherwise).
5. Report: the final code, how many attempts it took, the gate numbers (real TSs accepted, mirrored/noisy TSs caught),
   and share_of_goflow_guesses_it_flags. Be honest if the rule rarely fires: that means it adds little.
