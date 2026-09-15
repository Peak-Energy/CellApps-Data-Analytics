# lib/

Shared Python modules for logic used across more than one analysis category.
This is currently a scaffold — there are no real modules yet, since our
analyses run inside Voltaiq's analytics studio today. As we write code that's
reused in more than one place (or once we move analyses off Voltaiq), it
belongs here rather than being copy-pasted between scripts.

Organized by what the code *does*, not by data source, since most data
sources (Neware, Voltaiq, data logger, supplier files) will eventually need
the same kinds of processing/plotting regardless of where they came from:

- `io/` — loading and parsing raw data (Neware exports, Voltaiq pulls,
  data-logger files, the CaLT storage-time spreadsheet, supplier file
  formats, etc.)
- `processing/` — shared calculations (capacity/energy, DCIR, OCV extraction,
  round-trip efficiency, and similar)
- `plotting/` — shared chart styling and helpers, so figures look consistent
  across categories

If you write a function in a script and find yourself copying it into a
second script, that's the signal to move it here instead and import it in
both places.
