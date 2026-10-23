# Portable diagnostics

`diagnostics.json` is the small, machine-independent explanation artifact for one dataset. It is
written by preprocessing, candidate screening, standalone inference and `grnbench run`. Reports and
the interface can consume the same structure instead of reconstructing failures independently.

This is not a replacement for `manifest.json`: manifests retain exact paths, fingerprints,
dependencies and full provenance. Diagnostics are deliberately narrower and portable.

## Schema v1

Every document contains:

- `schema: "grnbench.diagnostics"` and `schema_version: 1`;
- the overall `outcome` enum (`unknown`, `ok`, `failed` or `partial`) and the
  original bounded `status` emitted by the stage;
- `prep`, `screen` and `models`, keyed by candidate/model and context;
- `usable_candidates` and `selected_candidate`.

Each candidate/context/model entry records:

- `ok` or `failed`;
- a stable failure category and message of at most 300 characters;
- applied choices; for preprocessing, `strategy` is the effective strategy and the optional
  `requested_strategy` preserves the original request (`auto`, for example);
- cells and genes before and after the relevant stage, when known;
- compact stage records for preprocessing.

The reader requires the complete schema envelope and accepts only the choice fields emitted by
the writer. Threshold choices are limited to the named QC thresholds and all other choice values
must be scalar. Missing schema fields, unknown choice fields and unsanitised messages are errors;
they are never passed through to a report.

Unknown dimensions are `null`; they are never guessed. A failure remains visible even when another
candidate or model succeeds.

## Privacy boundary

The strict writer recursively refuses:

- absolute paths and fields that contain paths;
- fields reserved for usernames, hostnames, email-like identity tokens, student IDs or barcodes;
- the current machine username and hostname when they occur as values;
- expression matrices or values;
- CLI commands and machine-specific output locations.

Exception messages are whitespace-normalized, known input/output paths, POSIX/Windows/UNC/file-URI
paths and local identity tokens are redacted, and the result is truncated to 300 characters before
validation. Candidate, context and model names are data, not schema field names, but their values
pass through the same path/identity validation. Raw data and detailed provenance remain in their
existing controlled artifacts; they are never copied into `diagnostics.json`.

Runtime integration is best effort: validation or serialization trouble emits a warning and never
turns an otherwise successful preprocessing, screening or inference run into a failure. The strict
writer remains available to tests and consumers that want an exception. JSON uses UTF-8, LF line
endings and no non-standard `NaN`/`Infinity` literals.

A directory containing only a valid schema-v1 diagnostic with `outcome: failed` may be reused for
the retry recommended by the CLI. Any other existing artifact still makes output reuse fail before
the command writes anything.

## Compatibility

Consumers must check both `schema` and `schema_version`. New optional entry fields are additive
within schema v1, so consumers must ignore fields they do not recognise. A future incompatible
shape increments the integer version; existing v1 documents remain readable. Registered experiment
outputs are not rewritten retroactively.

`grnbench.diagnostics.read_diagnostics` is the canonical reader. A missing file means a historical
run and is rendered explicitly as unavailable. An existing malformed, unsupported or privacy-unsafe
file raises an error; consumers never fall back to guessing its content from unrelated manifests.

Markdown, HTML and LaTeX reporting use the same normalized diagnostic rows. They show the overall
outcome, selected and usable candidates, and each preparation, screening and model outcome with its
context, choices, failure category/message, and before/after dimensions. `grnbench eval` accepts
`--diagnostics diagnostics.json` for a standalone evaluation; `grnbench run` refreshes its static
evaluation reports from the final root diagnostic automatically. Older reports remain generatable
and state why the diagnostic section is unavailable. The refreshed evaluation provenance
fingerprints that final file itself (not the surrounding run directory), and the one-shot command
prints the final report rather than its temporary pre-evaluation diagnostic.
