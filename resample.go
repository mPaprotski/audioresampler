// Package resample provides audio resampling functionality in pure Go.
// It converts audio signals from one sample rate to another using
// polyphase FIR filters with sinc interpolation.
package resample

import (
	"fmt"
	"math"
)

const (
	// defaultQuality controls the number of zero-crossings in the sinc filter.
	// Higher values = better quality but slower.
	defaultQuality = 64

	// defaultCutoffFactor is the fraction of Nyquist used as the low-pass cutoff.
	defaultCutoffFactor = 0.95
)

// Resampler holds precomputed filter coefficients for a specific
// input/output sample rate pair.
type Resampler struct {
	inRate  int
	outRate int
	// up/down are the reduced ratio numerator/denominator
	up      int
	down    int
	filters [][]float64 // polyphase filter bank [phase][tap]
	numTaps int
}

// New creates a Resampler for converting from inRate to outRate.
func New(inRate, outRate int) (*Resampler, error) {
	if inRate <= 0 || outRate <= 0 {
		return nil, fmt.Errorf("sample rates must be positive, got %d and %d", inRate, outRate)
	}
	g := gcd(inRate, outRate)
	up := outRate / g
	down := inRate / g

	r := &Resampler{
		inRate:  inRate,
		outRate: outRate,
		up:      up,
		down:    down,
	}
	r.buildFilterBank()
	return r, nil
}

// Resample converts audio from inRate to outRate.
// It is a convenience function that creates a temporary Resampler.
func Resample(source []int16, sampleRate int, targetSampleRate int) ([]int16, error) {
	if sampleRate == targetSampleRate {
		out := make([]int16, len(source))
		copy(out, source)
		return out, nil
	}
	r, err := New(sampleRate, targetSampleRate)
	if err != nil {
		return nil, err
	}
	return r.Process(source)
}

// Process resamples the given PCM samples.
func (r *Resampler) Process(input []int16) ([]int16, error) {
	if len(input) == 0 {
		return []int16{}, nil
	}

	// Convert int16 → float64 for processing
	in := make([]float64, len(input))
	for i, v := range input {
		in[i] = float64(v)
	}

	out := r.resampleFloat(in)

	// Convert float64 → int16 with clamping
	result := make([]int16, len(out))
	for i, v := range out {
		v = math.Round(v)
		if v > math.MaxInt16 {
			v = math.MaxInt16
		} else if v < math.MinInt16 {
			v = math.MinInt16
		}
		result[i] = int16(v)
	}
	return result, nil
}

// resampleFloat performs polyphase resampling on float64 samples.
// Algorithm:
//  1. Upsample by factor `up` (insert up-1 zeros between samples)
//  2. Low-pass filter to remove imaging/aliasing
//  3. Downsample by factor `down` (keep every down-th sample)
//
// Steps 1-3 are combined efficiently using a polyphase decomposition
// so we never actually materialise the upsampled signal.
func (r *Resampler) resampleFloat(input []float64) []float64 {
	numTaps := r.numTaps
	halfTaps := numTaps / 2

	// Pad input with zeros on both sides so the filter is causal
	padded := make([]float64, halfTaps+len(input)+halfTaps)
	copy(padded[halfTaps:], input)

	// Compute output length: ceil(len(input) * up / down)
	outLen := (len(input)*r.up + r.down - 1) / r.down
	output := make([]float64, outLen)

	for i := 0; i < outLen; i++ {
		// Output sample i corresponds to input position i*down/up
		// Phase within the polyphase filter bank
		phase := (i * r.down) % r.up
		// Input sample index (centre tap)
		inputIdx := (i*r.down)/r.up + halfTaps

		filter := r.filters[phase]
		var acc float64
		for j, coef := range filter {
			srcIdx := inputIdx - halfTaps + j
			if srcIdx >= 0 && srcIdx < len(padded) {
				acc += padded[srcIdx] * coef
			}
		}
		// Scale by upsampling factor (energy compensation)
		output[i] = acc * float64(r.up)
	}
	return output
}

// buildFilterBank precomputes a polyphase FIR filter bank.
// The prototype filter is a windowed-sinc low-pass filter.
func (r *Resampler) buildFilterBank() {
	quality := defaultQuality

	// Cutoff is min(1/up, 1/down) * cutoffFactor, normalised to [0,1] where 1 = Nyquist
	cutoff := defaultCutoffFactor / float64(max(r.up, r.down))

	// Total taps in prototype filter (must be odd for linear phase)
	totalTaps := 2*quality*max(r.up, r.down) + 1
	r.numTaps = totalTaps

	// Build prototype windowed-sinc filter
	proto := make([]float64, totalTaps)
	mid := (totalTaps - 1) / 2
	for i := 0; i < totalTaps; i++ {
		n := float64(i - mid)
		// Sinc
		var sinc float64
		if n == 0 {
			sinc = 2 * cutoff
		} else {
			sinc = math.Sin(2*math.Pi*cutoff*n) / (math.Pi * n)
		}
		// Kaiser window (β=6) for good stopband attenuation
		proto[i] = sinc * kaiserWindow(i, totalTaps, 6.0)
	}

	// Decompose into r.up polyphase sub-filters
	r.filters = make([][]float64, r.up)
	for phase := 0; phase < r.up; phase++ {
		// Each sub-filter picks every r.up-th coefficient starting at `phase`
		subLen := 0
		for k := phase; k < totalTaps; k += r.up {
			subLen++
		}
		sub := make([]float64, subLen)
		idx := 0
		for k := phase; k < totalTaps; k += r.up {
			sub[idx] = proto[k]
			idx++
		}
		r.filters[phase] = sub
	}
	// Update numTaps to sub-filter length (used for padding)
	r.numTaps = len(r.filters[0])
}

// kaiserWindow computes the Kaiser window value at position n
// for a window of length M with shape parameter beta.
func kaiserWindow(n, m int, beta float64) float64 {
	x := 2.0*float64(n)/float64(m-1) - 1.0
	return besselI0(beta*math.Sqrt(math.Max(0, 1-x*x))) / besselI0(beta)
}

// besselI0 computes the modified Bessel function of the first kind, order 0.
// Uses a polynomial approximation accurate to ~1e-9.
func besselI0(x float64) float64 {
	ax := math.Abs(x)
	if ax < 3.75 {
		y := (x / 3.75) * (x / 3.75)
		return 1.0 + y*(3.5156229+y*(3.0899424+y*(1.2067492+
			y*(0.2659732+y*(0.0360768+y*0.0045813)))))
	}
	y := 3.75 / ax
	return (math.Exp(ax) / math.Sqrt(ax)) * (0.39894228 + y*(0.01328592+
		y*(0.00225319+y*(-0.00157565+y*(0.00916281+y*(-0.02057706+
			y*(0.02635537+y*(-0.01647633+y*0.00392377))))))))
}

// gcd returns the greatest common divisor of a and b.
func gcd(a, b int) int {
	for b != 0 {
		a, b = b, a%b
	}
	return a
}

// max returns the larger of a and b.
func max(a, b int) int {
	if a > b {
		return a
	}
	return b
}
