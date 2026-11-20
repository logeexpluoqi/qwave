/**
 * Author: luoqi
 * Created Date: 2025-12-23 16:08:59
 * Last Modified: 2026-11-21 00:37:39
 * Modified By: luoqi at <**@****>
 * Copyright (c) 2026 <*****>
 * Description:
 */

#ifndef _QWAVE_H_
#define _QWAVE_H_

#include <stdint.h>

#ifndef fp_t
typedef float fp_t;
#endif

#ifndef NAN
#define NAN (0.0f / 0.0f)
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
    QWaveType type;      /**< Type of the waveform. */
    fp_t fs;             /**< Sampling frequency in Hz. Must be greater than 5 times the wave frequency. */
    fp_t frq;            /**< Wave frequency in Hz. */
    fp_t period;         /**< Period of the waveform in seconds. */
    fp_t half_period;    /**< Half of the period in seconds. */
    fp_t t;              /**< Current time or phase of the waveform. */
    fp_t ts;             /**< Time step or sampling interval. */
    fp_t bias;           /**< Bias or offset of the waveform. */
    fp_t amp;            /**< Amplitude of the waveform. */
    fp_t output;         /**< Current output value of the waveform. */
    uint32_t seed;       /**< Seed value for the pseudo-random number generator. */
    uint32_t prng_state; /**< State of the pseudo-random number generator. */
} QWaveCore;

/**
 * @brief Initialize a QWaveCore structure with the given parameters, and set the initial output value.
 * @param core: pointer to the QWaveCore structure to be initialized.
 * @param type: type of the waveform (QWaveType).
 * @param fs: sampling frequency in Hz.
 * @param frq: wave frequency in Hz.
 * @param bias: bias or offset of the waveform.
 * @param seed: seed value for the pseudo-random number generator.
 * @return 0 on success, non-zero on failure.
 */
int qwave_init(QWaveCore *core, QWaveType type, fp_t fs, fp_t frq, fp_t bias, uint32_t seed);

/**
 * @brief Generate the next output value of the waveform based on the current state.
 * @param core: pointer to the QWaveCore structure.
 * @return The next output value of the waveform.
 * @note This function updates the internal state of the waveform generator and returns the next output value,
 *       sampling frequency is depending on parameter fs.
 */
fp_t qwave_output(QWaveCore *core);

/**
 * @brief Set the bias or offset of the waveform.
 * @param core: pointer to the QWaveCore structure.
 * @param bias: new bias or offset value.
 * @return 0 on success, non-zero on failure.
 */
int qwave_bias_set(QWaveCore *core, fp_t bias);

/**
 * @brief Set the sampling frequency of the waveform.
 * @param core: pointer to the QWaveCore structure.
 * @param fs: new sampling frequency in Hz.
 * @return 0 on success, non-zero on failure.
 */
int qwave_fs_set(QWaveCore *core, fp_t fs);

/**
 * @brief Set the wave frequency of the waveform.
 * @param core: pointer to the QWaveCore structure.
 * @param frq: new wave frequency in Hz.
 * @return 0 on success, non-zero on failure.
 */
int qwave_frq_set(QWaveCore *core, fp_t frq);

/**
 * @brief Set the amplitude of the waveform.
 * @param core: pointer to the QWaveCore structure.
 * @param amp: new amplitude value.
 * @return 0 on success, non-zero on failure.
 */
int qwave_amp_set(QWaveCore *core, fp_t amp);

/**
 * @brief Set the type of the waveform.
 * @param core: pointer to the QWaveCore structure.
 * @param type: new type of the waveform (QWaveType).
 * @return 0 on success, non-zero on failure.
 */
int qwave_signal_set(QWaveCore *core, QWaveType type);

/**
 * @brief Reset the QWaveCore structure, resetting it to its initial state.
 * @param core: pointer to the QWaveCore structure.
 * @return 0 on success, non-zero on failure.
 */
int qwave_reset(QWaveCore *core);

#endif
