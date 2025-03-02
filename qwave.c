/**
 * @ Author: luoqi
 * @ Create Time: 2024-11-08 17:16
 * @ Modified by: luoqi
 * @ Modified time: 2025-03-02 20:08
 * @ Description:
 */

#include "qwave.h"

static const qfp_t _fast_sin_table[91] = {
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

static inline qfp_t _rpm2deg(qfp_t rpm)
{
    return 6 * rpm;
}

static inline qfp_t _deg2rpm(qfp_t deg)
{
    return deg / 6;
}

static inline qfp_t _fmodf(qfp_t x, qfp_t y)
{
    if(y == 0) {
        return NAN;
    }
    return (x - ((int)(x / y)) * y);
}

static qfp_t _fsin(qfp_t deg)
{
    deg = _fmodf(deg, 360);

    int sign = 1;
    // Use the periodic symmetry of sin
    if(deg > 180) {
        deg -= 180;
        sign = -1;
    }
    if(deg > 90) {
        deg = 180 - deg;
    }

    // Calculate the lookup table index and interpolation factor
    int index = (int)deg;  // 1° resolution
    if(index >= 90) {
        return sign * _fast_sin_table[90];
    }
    qfp_t fraction = deg - index;

    qfp_t sin_val = _fast_sin_table[index] * (1 - fraction) + _fast_sin_table[index + 1] * fraction;
    return sign * sin_val;
}

static inline qfp_t _fcos(qfp_t x)
{
    return _fsin(x + 90);
}

static inline qfp_t _gen_sin(QWaveGen *gen)
{
    qfp_t x = (gen->t / gen->period) * 360;
    gen->output = _fsin(gen->frq * x) + gen->bias;
    return gen->output;
}

static inline qfp_t _gen_tri(QWaveGen *gen)
{
    if(!gen) {
        return 0;
    }

    qfp_t norm = gen->t / gen->period;

    if(norm < 0.25) {
        gen->output = 4 * norm;
    } else if(norm < 0.75) {
        gen->output = 2 - 4 * norm;
    } else {
        gen->output = 4 * norm - 4;
    }
    gen->output += gen->bias;
    return gen->output;
}

static inline qfp_t _gen_saw(QWaveGen *gen)
{
    if(!gen) {
        return 0;
    }
    gen->output = gen->t / gen->period;
    return gen->output;
}

static inline qfp_t _gen_antsaw(QWaveGen *gen)
{
    if(!gen) {
        return 0;
    }
    gen->output = - (gen->t / gen->period);
    return gen->output;
}

static inline qfp_t _gen_sqr(QWaveGen *gen)
{
    if(!gen) {
        return 0;
    }
    // Square wave: first half period: +1; second half: -1.
    if(gen->t < gen->half_period) {
        gen->output = 1;
    } else {
        gen->output = -1;
    }
    gen->output += gen->bias;
    return gen->output;
}

static inline qfp_t _gen_noise(QWaveGen *gen)
{
    if(!gen) {
        return 0;
    }
    uint32_t x = gen->prng_state;
    x ^= x << 13;
    x ^= x >> 17;
    x ^= x << 5;
    gen->prng_state = x;
    // Map x (0 ~ 0xffffffffu) to [-1, 1]
    gen->output = ((qfp_t)x / 0xffffffffu) * 2 - 1 + gen->bias;
    return gen->output;
}

int qwave_tick_init(QWaveGen *gen, QWaveType type, qfp_t fs, qfp_t frq, qfp_t bias, uint32_t seed)
{
    if(!gen || fs <= 0 || frq <= 0) {
        return -1;
    }

    gen->type = type;
    gen->bias = bias;
    gen->fs = fs;
    gen->frq = frq;
    gen->period = 1.0 / frq;
    gen->half_period = gen->period * 0.5;
    gen->ts = 1.0 / fs;
    gen->output = 0;
    gen->amp = 0.0;
    gen->t = 0;
    gen->seed = seed;
    gen->prng_state = (seed == 0) ? 2463534242UL : seed;

    return 0;
}

int qwave_time_init(QWaveGen *gen, QWaveType type, qfp_t frq, qfp_t bias, uint32_t seed)
{
    return qwave_tick_init(gen, type, 10 * frq, frq, bias, seed);
}

int qwave_signal_set(QWaveGen *gen, QWaveType type)
{
    if(!gen) {
        return -1;
    }
    gen->type = type;
    gen->t = 0;
    gen->output = 0;
    return 0;
}

static inline qfp_t _qwave_out(QWaveGen *gen)
{
    switch(gen->type) {
    case QWAVE_TYPE_SINE:
        return _gen_sin(gen);
    case QWAVE_TYPE_TRIANGLE:
        return _gen_tri(gen);
    case QWAVE_TYPE_SAWTOOTH:
        return _gen_saw(gen);
    case QWAVE_TYPE_ANTSAWTOOTH:
        return _gen_antsaw(gen);
    case QWAVE_TYPE_SQUARE:
        return _gen_sqr(gen);
    case QWAVE_TYPE_NOISE:
        return _gen_noise(gen);
    default:
        return 0;
    }
}

qfp_t qwave_tick_signal_output(QWaveGen *gen)
{
    if(!gen) {
        return 0;
    }

    qfp_t out = _qwave_out(gen);

    gen->t += gen->ts;
    if(gen->t >= gen->period) {
        gen->t -= gen->period;
    }
    return out;
}

qfp_t qwave_time_signal_output(QWaveGen *gen, qfp_t dms)
{
    if(!gen) {
        return 0;
    }
    qfp_t out = _qwave_out(gen) * gen->amp;
    gen->t += dms * 1e-3;
    if(gen->t >= gen->period) {
        gen->t -= gen->period;
    }
    return out;
}

int qwave_bias_set(QWaveGen *gen, qfp_t bias)
{
    if(!gen) {
        return -1;
    }
    gen->bias = bias;
    return 0;
}

int qwave_fs_set(QWaveGen *gen, qfp_t fs)
{
    if(!gen || fs <= 0) {
        return -1;
    }
    gen->fs = fs;
    gen->ts = 1.0 / fs;
    return 0;
}

int qwave_frq_set(QWaveGen *gen, qfp_t frq)
{
    if(!gen || frq <= 0) {
        return -1;
    }
    gen->frq = frq;
    gen->period = 1.0 / frq;
    gen->half_period = gen->period * 0.5;
    return 0;
}

int qwave_amp_set(QWaveGen *gen, qfp_t amp)
{
    if(!gen || amp <= 0) {
        return -1;
    }
    gen->amp = amp;
    return 0;
}

int qwave_clr(QWaveGen *gen)
{
    if(!gen) {
        return -1;
    }
    gen->t = 0;
    gen->output = 0;
    return 0;
}
