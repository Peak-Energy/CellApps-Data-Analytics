# Battery Data Analytics

Central repository for Peak Energy's battery data analysis code. This repo is
organized by data category rather than by tool, so that as our workflow
evolves (e.g. moving analyses out of Voltaiq over time) the folder structure
doesn't need to change.

## What lives here (and what doesn't)

**Raw data is never committed to this repo.** Test data lives in Voltaiq,
on SharePoint, or on local data loggers, depending on the category — see each
category's own README for specifics. This repo holds only code:

- **Voltaiq analytics studio templates** (`voltaiq_scripts/` folders): most of
  our current analysis runs inside Voltaiq's analytics studio, not locally.
  The scripts checked in here are template/reference copies of what's running
  in Voltaiq, kept for version history and so the whole team has a single
  place to find them. They are not necessarily meant to be run standalone.
- **Local scripts** (`scripts/` folders): code that actually runs locally
  against data files (e.g. supplier data received by email, or data-logger
  output for Entropy testing).
- **Shared library code** (`lib/`): reusable modules for parsing, processing,
  and plotting. This is currently a scaffold with no real modules yet — as we
  write code that's used across more than one analysis, it belongs here. See
  `lib/README.md`.

## Repository layout

```
baseline/           # CyLT, CaLT, DCIR/OCV — Neware-run, tracked in Voltaiq
detailed/           # HPPC, FastGITT (Voltaiq), Entropy (local data logger)
thermal/            # Energy Remaining, ARC, IBC
mechanical/         # Swell Force Test
supplier/           # Supplier-provided data: Veken, MCM, EVE, EA
lib/                # Shared code: io/, processing/, plotting/ (scaffold for now)
```

Each leaf folder (e.g. `baseline/cylt/`) has its own README describing what
the analysis covers, where its source data lives, and any quirks specific to
that category.

## Getting started

Clone the repo to get a local copy of the current scripts:

```
git clone <repo-url>
```

There's no installable package or dependency file yet since nothing here
runs locally besides the `supplier/` and `detailed/entropy/` scripts. If you
add a script with real Python dependencies, add a `requirements.txt` (or
`environment.yml`) alongside it, or at the repo root if it's shared.

## Contributing

Branch off `main` for any change (`git checkout -b your-branch-name`), and
open a pull request when it's ready for another team member to look at
before merging. See individual category READMEs before adding or editing a
script in that folder.
