# resample

Pure-Go audio resampling library — no C dependencies, no cgo.

Designed for use in speech recognition, speech synthesis, and voice biometrics pipelines.

## Features

| Feature | Detail |
|---------|--------|
| Integer upsampling | 8 000 → 16 000 Hz (×2), 8 000 → 48 000 Hz (×6), … |
| Integer downsampling | 48 000 → 8 000 Hz (÷6), 44 100 → 22 050 Hz (÷2), … |
| Fractional resampling | 11 025 → 8 000 Hz (×160/441), 11 025 → 16 000 Hz, … |
| Anti-aliasing | Windowed-sinc FIR with Kaiser window (β = 6) |
| No dependencies | Pure Go, zero cgo |

## Quick start

```go
import "github.com/example/resample"

// One-shot convenience function
output, err := resample.Resample(input, 8000, 16000)

// Reusable resampler (filter built once, then reused)
r, err := resample.New(8000, 16000)
output, err = r.Process(chunk1)
output2, err = r.Process(chunk2)
```

### Signature

```go
func Resample(source []int16, sampleRate int, targetSampleRate int) ([]int16, error)
```

## Algorithm

The library uses **polyphase FIR resampling**:

1. Represent the rational conversion factor as `up/down` (reduced by GCD).
2. Build a prototype low-pass FIR filter (windowed-sinc, Kaiser window β = 6)
   with cutoff at `min(1/up, 1/down) × 0.95 × Nyquist`.
3. Decompose the prototype into `up` polyphase sub-filters.
4. For each output sample: select the correct sub-filter phase and compute
   the dot product with the relevant input samples.

This avoids ever materialising the ×up upsampled signal in memory and runs in
`O(len(input) × numTaps / down)` time.

## Running tests

```bash
go test ./...
go test -bench=. -benchmem ./...
```

## Common conversions

| From → To | Ratio |
|-----------|-------|
| 8 000 → 16 000 | ×2 |
| 16 000 → 8 000 | ÷2 |
| 48 000 → 8 000 | ÷6 |
| 48 000 → 16 000 | ÷3 |
| 11 025 → 8 000 | ×160/441 |
| 11 025 → 16 000 | ×320/441 |
| 44 100 → 16 000 | ×160/441 |