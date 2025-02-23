/**
 * @ Author: luoqi
 * @ Create Time: 2024-11-08 17:16
 * @ Modified by: luoqi
 * @ Modified time: 2025-02-23 23:51
 * @ Description:
 */

#ifndef _QWAVE_H
#define _QWAVE_H

#include <stdint.h>

#ifndef qfp_t
typedef float qfp_t;
#endif

typedef enum {
    QWAVE_TYPE_SINE = 0,
    QWAVE_TYPE_TRIANGLE,
    QWAVE_TYPE_SAWTOOTH,
    QWAVE_TYPE_ANTSAWTOOTH,
    QWAVE_TYPE_NOISE,
    QWAVE_TYPE_SQUARE
} QWaveType;

typedef struct {
    QWaveType type;       /**< Type of the waveform. */
    qfp_t fs;             /**< Sampling frequency in Hz. Must be greater than 5 times the wave frequency. */
    qfp_t frq;            /**< Wave frequency in Hz. */
    qfp_t period;         /**< Period of the waveform in seconds. */
    qfp_t half_period;    /**< Half of the period in seconds. */
    qfp_t t;              /**< Current time or phase of the waveform. */
    qfp_t ts;             /**< Time step or sampling interval. */
    qfp_t bias;           /**< Bias or offset of the waveform. */
    qfp_t output;         /**< Current output value of the waveform. */
    qfp_t timer;          /**< Timer value for internal use. */
    uint32_t seed;        /**< Seed value for the pseudo-random number generator. */
    uint32_t prng_state;  /**< State of the pseudo-random number generator. */
} QWaveGen;

int qwave_init(QWaveGen *gen, QWaveType type, qfp_t fs, qfp_t frq, qfp_t bias, uint32_t seed);

qfp_t qwave_generate(QWaveGen *gen);

int qwave_bias_set(QWaveGen *gen, qfp_t bias);

int qwave_fs_set(QWaveGen *gen, qfp_t fs);

int qwave_frq_set(QWaveGen *gen, qfp_t frq);

int qwave_signal_set(QWaveGen *gen, QWaveType type);

#endif
