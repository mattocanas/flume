# Contributing to Flume

Thanks for your interest in improving Flume! Contributions of all kinds are
welcome — bug reports, feature ideas, documentation, and code.

## Ground rules

- Be respectful and constructive.
- By contributing, you agree that your contributions are licensed under the
  project's [BSD 3-Clause License](./LICENSE).

## Reporting bugs / requesting features

Open an issue. For bugs, please include:

- What you did and what you expected vs. what happened.
- Your browser and OS.
- A small example file if the issue is data-specific (please ensure it contains
  no sensitive or unpublished data).

## Development setup

Requires [Node.js](https://nodejs.org/) 18 or newer.

```bash
git clone https://github.com/mattocanas/flume.git
cd flume
npm install
npm run dev
```

The app lives almost entirely in `src/flow_combined.jsx`; `src/fcs.js` handles
raw `.fcs` parsing. There are no third-party plotting libraries — all rendering
is hand-written Canvas/SVG.

## Pull requests

1. Fork and create a branch off `main`.
2. Keep changes focused; describe what and why in the PR.
3. Make sure `npm run build` succeeds before submitting.
4. Prefer matching the surrounding code style.

## Design principle

Flume is **100% client-side by design** — no accounts, no server, no data leaving
the user's machine. Please keep it that way: contributions must not add
network calls that transmit user data.
