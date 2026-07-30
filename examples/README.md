# Example data

Synthetic datasets for trying out Flume — **not real experimental data**. They're
generated to resemble a typical binding assay (viability × antibody signal) so
every plot type and stat has something meaningful to show.

| File | Events | Description |
|------|--------|-------------|
| `example_stained.csv` | 5,000 | Stained sample — a clear antibody-positive population on `Alexa Fluor 647-A`. |
| `example_isotype.csv` | 5,000 | Isotype control — mostly negative, for overlay/log₂ fold-change comparisons. |

**Channels:** `FSC-A`, `SSC-A`, `Zombie-Violet-A` (viability), `Alexa Fluor 647-A` (signal).

## Try it

1. Open the app and drag **both** CSVs onto the upload box.
2. **Histogram / Overlay / Ridge:** pick `Alexa Fluor 647-A`; drag the gate to
   split negative vs. positive.
3. **Quadrant:** set X = `Zombie-Violet-A`, Y = `Alexa Fluor 647-A` and drag the
   crosshair to gate live/dead × signal.
4. Compare the stained sample's gated gMFI and log₂ fold-change against the
   isotype control.
