# Outcome Templates

## 1. Run Completion Summary

## Objective
- what this run aimed to complete

## Result
- success / partial / failed

## Canonical evidence
- run directory
- `module_status.csv`
- `run_manifest.json`
- key output artifacts

## Main lesson
- the single most important carry-forward conclusion

## Next action
- the most useful next operational step

## 2. Failure Summary

## Failure surface
- first failing module or launch stage

## Observed error
- exact error text or concise paraphrase

## Root cause
- observed / inferred split if needed

## Recovery action
- what changed next

## Status
- fixed / not yet fixed / workaround only

## 3. Module Blocker Summary

## Module
- module name

## Required inputs
- what the module expected

## Blocker
- what failed first

## Fix applied
- code or execution change

## Proof of fix
- later run or artifact that proved the blocker moved

## 4. Article Claim Update Summary

## Claim family
- which paper figure / section / abstract claim

## Current support
- direct / partial / unsupported

## Our evidence
- concrete artifact paths

## Main gap
- what still prevents a stronger verdict

## 5. Cleanup Summary

## Preserved
- canonical runs / logs / inputs

## Deleted
- clearly superseded or partial artifacts

## Reasoning
- why deletion was safe

## Space impact
- approximate reclaimed size
