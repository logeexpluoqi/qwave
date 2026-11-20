/**
 * Author: luoqi
 * Created Date: 2025-12-23 16:08:59
 * Last Modified: 2026-11-21 00:37:39
 * Modified By: luoqi at <**@****>
 * Copyright (c) 2026 <*****>
 * Description:
 */

#include "qwave.h"

// clang-format off
static const fp_t fast_sin_table_[91] = {
    0.0,           0.017452406,   0.034899497,   0.052335956,   0.069756474,
    0.087155743,   0.104528463,   0.121869343,   0.139173101,   0.156434465,
    0.173648178,   0.190809,      0.207911,      0.224951,      0.241922,
    0.258819,      0.275637,      0.292372,      0.309017,      0.325568,
    0.342020,      0.358368,      0.374607,      0.390731,      0.406737,
    0.422618,      0.438371,      0.453990,      0.469472,      0.48481,
    0.5,           0.515038,      0.529919,      0.544639,      0.559193,
    0.573576,      0.587785,      0.601815,      0.615661,      0.62932,
    0.642788,      0.656059,      0.669131,      0.681998,      0.694658,
    0.707107,      0.71934,       0.731354,      0.743145,      0.75471,
    0.766044,      0.777146,      0.788011,      0.798636,      0.809017,
    0.819152,      0.829038,      0.838671,      0.848048,      0.857167,
    0.866025,      0.87462,       0.882948,      0.891007,      0.898794,
    0.906308,      0.913545,      0.920505,      0.927184,      0.93358,
    0.939693,      0.945519,      0.951057,      0.956305,      0.961262,
    0.965926,      0.970296,      0.97437,       0.978148,      0.981627,
    0.984808,      0.987688,      0.990268,      0.992546,      0.994522,
    0.996195,      0.997564,      0.99863,       0.999391,      0.999848,
    1.0
};
// clang-format on

static inline fp_t fmodf_(fp_t x, fp_t y)
{
    if(y == 0) {
        return NAN;
    }

    fp_t abs_y = y < 0 ? -y : y;
    fp_t sign_x = x < 0 ? -1.0f : 1.0f;
    fp_t abs_x = x < 0 ? -x : x;

    fp_t div = abs_x / abs_y;
    fp_t int_part = div >= 0 ? (fp_t)((int)div) : -(fp_t)((int)-div);

    fp_t result = abs_x - (int_part * abs_y);

    return sign_x * result;
}

static fp_t sin_(fp_t deg)
{
    deg = fmodf_(deg, 360.0f);

    int sign = 1;
    // Use the periodic symmetry of sin
    if(deg > 180.0f) {
        deg -= 180.0f;
        sign = -1;
    }
    if(deg > 90.0f) {
        deg = 180.0f - deg;
    }

    // Calculate the lookup table index and interpolation factor
    int index = (int)deg; // 1° resolution
    if(index >= 90) {
        return sign * fast_sin_table_[90];
    }
    fp_t fraction = deg - index;

    fp_t sin_val = fast_sin_table_[index] * (1.0f - fraction) + fast_sin_table_[index + 1] * fraction;
    return sign * sin_val;
}

/* All gen_*_ functions return the raw waveform in normalized units.
 * Amplitude scaling and bias are applied centrally in qwave_output(),
 * so that amp never scales the DC bias. */

static inline fp_t gen_sin_(QWaveCore *core)
{
    return sin_((core->t * 360.0f) * core->frq);
}

static inline fp_t gen_tri_(QWaveCore *core)
{
    fp_t norm = core->t * core->frq;

    if(norm < 0.25f) {
        return 4.0f * norm;
    } else if(norm < 0.75f) {
        return 2.0f - 4.0f * norm;
    }
    return 4.0f * norm - 4.0f;
}

static inline fp_t gen_saw_(QWaveCore *core)
{
    return core->t * core->frq;
}

static inline fp_t gen_antsaw_(QWaveCore *core)
{
    return -(core->t * core->frq);
}

static inline fp_t gen_sqr_(QWaveCore *core)
{
    // Square wave: first half period: +1; second half: -1.
    return (core->t < core->half_period) ? 1.0f : -1.0f;
}

static inline fp_t gen_noise_(QWaveCore *core)
{
    uint32_t x = core->prng_state;
    x ^= x << 13;
    x ^= x >> 17;
    x ^= x << 5;
    core->prng_state = x;
    // Map x (0 ~ 0xffffffffu) to [-1, 1]
    return ((fp_t)x / (fp_t)0xffffffffu) * 2.0f - 1.0f;
}

int qwave_init(QWaveCore *core, QWaveType type, fp_t fs, fp_t frq, fp_t bias, uint32_t seed)
{
    if(!core || fs <= 0 || frq <= 0 || type < QWAVE_TYPE_SINE || type > QWAVE_TYPE_SQUARE) {
        return -1;
    }

    core->type = type;
    core->bias = bias;
    core->fs = fs;
    core->frq = frq;
    core->period = 1.0f / frq;
    core->half_period = core->period * 0.5f;
    core->ts = 1.0f / fs;
    core->output = 0;
    core->amp = 1.0f;
    core->t = 0;
    core->seed = seed;
    core->prng_state = (seed == 0) ? 2463534242UL : seed;

    return 0;
}

int qwave_signal_set(QWaveCore *core, QWaveType type)
{
    if(!core || type < QWAVE_TYPE_SINE || type > QWAVE_TYPE_SQUARE) {
        return -1;
    }
    core->type = type;
    core->t = 0;
    core->output = 0;
    return 0;
}

static inline fp_t qwave_out_(QWaveCore *core)
{
    switch(core->type) {
    case QWAVE_TYPE_SINE:
        return gen_sin_(core);
    case QWAVE_TYPE_TRIANGLE:
        return gen_tri_(core);
    case QWAVE_TYPE_SAWTOOTH:
        return gen_saw_(core);
    case QWAVE_TYPE_ANTSAWTOOTH:
        return gen_antsaw_(core);
    case QWAVE_TYPE_SQUARE:
        return gen_sqr_(core);
    case QWAVE_TYPE_NOISE:
        return gen_noise_(core);
    default:
        return 0;
    }
}

fp_t qwave_output(QWaveCore *core)
{
    if(!core) {
        return 0;
    }

    /* amp scales only the raw waveform; bias is a pure DC offset */
    core->output = qwave_out_(core) * core->amp + core->bias;

    core->t += core->ts;
    core->t = fmodf_(core->t, core->period);
    return core->output;
}

int qwave_bias_set(QWaveCore *core, fp_t bias)
{
    if(!core) {
        return -1;
    }
    core->bias = bias;
    return 0;
}

int qwave_fs_set(QWaveCore *core, fp_t fs)
{
    if(!core || fs <= 0) {
        return -1;
    }
    core->fs = fs;
    core->ts = 1.0f / fs;
    return 0;
}

int qwave_frq_set(QWaveCore *core, fp_t frq)
{
    if(!core || frq <= 0) {
        return -1;
    }
    core->frq = frq;
    core->period = 1.0f / frq;
    core->half_period = core->period * 0.5f;
    /* wrap the phase into the new period to avoid a one-sample phase jump */
    core->t = fmodf_(core->t, core->period);
    return 0;
}

int qwave_amp_set(QWaveCore *core, fp_t amp)
{
    if(!core || amp < 0) {
        return -1;
    }
    core->amp = amp;
    return 0;
}

int qwave_reset(QWaveCore *core)
{
    if(!core) {
        return -1;
    }
    core->t = 0;
    core->output = 0;
    /* restore the PRNG so noise sequences are reproducible after reset */
    core->prng_state = (core->seed == 0) ? 2463534242UL : core->seed;
    return 0;
}
