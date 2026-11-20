# qwave

qwave 是一个轻量级 C 波形发生器库。它以固定采样率逐点生成周期信号或伪随机噪声，适用于嵌入式程序、信号处理原型和仿真代码。库仅依赖标准 C 头文件，不需要额外的数学库。

## 特性

- 正弦波、三角波、锯齿波、反向锯齿波、方波和伪随机噪声。
- 振幅和直流偏置分开控制，振幅不会缩放偏置。
- 可在运行时修改波形类型、采样率、频率、振幅和偏置。
- 噪声使用可复现的 xorshift32 伪随机序列；重置后可从相同种子重新生成。
- 每个发生器的状态都保存在 `QWaveCore` 实例中，可同时使用多个独立实例。

## 文件

- `qwave.h`：类型定义、波形枚举和公开 API。
- `qwave.c`：波形计算及状态更新实现。
- `sim.m`：MATLAB/Octave 信号生成仿真脚本，不参与 C 库构建。

## 构建

将 `qwave.c` 与使用它的 C 源文件一起编译，并将项目目录加入头文件搜索路径。例如：

```sh
cc -std=c99 -O2 -Wall -Wextra -I. main.c qwave.c -o main
```

无需链接 `libm`。如需将库集成到现有工程，把 `qwave.c` 编译进工程并在调用处包含 `qwave.h` 即可。

## 快速开始

```c
#include <stdio.h>
#include "qwave.h"

int main(void)
{
    QWaveCore core;

    /* 1 kHz 采样，100 Hz 正弦波，零偏置；seed 对非噪声波形无影响。 */
    if(qwave_init(&core, QWAVE_TYPE_SINE, 1000.0f, 100.0f, 0.0f, 1) != 0) {
        return 1;
    }

    qwave_amp_set(&core, 2.0f);
    qwave_bias_set(&core, 0.5f);

    for(int i = 0; i < 10; ++i) {
        printf("%f\n", qwave_output(&core));
    }

    return 0;
}
```

每次调用 `qwave_output()` 会生成并返回一个样本，同时推进发生器状态。上例的采样间隔为 `1 / fs = 1 ms`。首个样本在相位 `t = 0` 处生成。

## 设计与实现

### 统一输出模型

各波形生成函数先产生归一化原始值 `raw`，再由 `qwave_output()` 统一计算：

```text
output = raw * amp + bias
```

因此 `amp` 只缩放波形本身，`bias` 始终作为独立的直流偏置。初始化时 `amp` 默认为 `1`，偏置取 `qwave_init()` 传入的值。比如 `amp = 2`、`bias = 0.5` 时，正弦输出范围为 `[-1.5, 2.5]`。

### 相位和采样

发生器用 `t`（秒）表示当前相位时间，用 `ts = 1 / fs` 表示采样间隔。每次输出按以下顺序工作：

1. 根据当前 `t` 和波形类型计算样本。
2. 应用振幅和偏置并保存到 `core->output`。
3. 将 `t` 增加 `ts`，并按周期 `period = 1 / frq` 回绕。

周期波形以该周期重复。采样率用于控制相邻样本的时间间隔；初始化和采样率设置要求 `fs > 0`，但库不会检查采样率与波形频率之间的比例。头文件建议 `fs` 大于 `5 * frq`；实际使用时应根据目标信号和下游处理要求选择采样率，避免混叠。

### 各波形的原始值

| 枚举 | 原始波形 | 说明 |
| --- | --- | --- |
| `QWAVE_TYPE_SINE` | `sin(2π · frq · t)` | 正弦值通过 1° 分辨率查找表和线性插值得到，不调用系统三角函数。 |
| `QWAVE_TYPE_TRIANGLE` | 周期内从 `0` 上升至 `+1`，再下降至 `-1`，最后回到 `0` | 起始相位为零，四分之一周期到达正峰值，四分之三周期到达负峰值。 |
| `QWAVE_TYPE_SAWTOOTH` | `t / period` | 每周期从 `0` 线性上升到接近 `+1`，随后回到 `0`。 |
| `QWAVE_TYPE_ANTSAWTOOTH` | `-t / period` | 每周期从 `0` 线性下降到接近 `-1`，随后回到 `0`。 |
| `QWAVE_TYPE_SQUARE` | 前半周期 `+1`，后半周期 `-1` | 在周期中点切换；新周期开始时回到 `+1`。 |
| `QWAVE_TYPE_NOISE` | `[-1, +1]` 上的伪随机值 | 每次调用推进一次 xorshift32 状态；噪声值不依赖 `t`，但采样仍按 `fs` 调用。 |

正弦波使用固定表降低计算开销，线性插值处理整数角度间的输入。锯齿类和三角波基于归一化周期位置计算；方波依据半周期判断。上述波形均先在 `[-1, 1]` 或 `[0, 1)` 等归一化范围内生成，再应用统一输出模型。

## API 参考

包含头文件：

```c
#include "qwave.h"
```

### 类型

```c
typedef float fp_t;
```

`fp_t` 是库内部的浮点类型，默认定义为 `float`。`QWaveCore` 是公开状态结构体，包含波形类型、采样率、频率、周期、当前时间、采样间隔、偏置、振幅、最近输出以及噪声随机状态。通常应通过库函数操作状态，而不是直接改写结构体字段。

波形枚举值：

```c
QWAVE_TYPE_SINE          /* 正弦波 */
QWAVE_TYPE_TRIANGLE      /* 三角波 */
QWAVE_TYPE_SAWTOOTH      /* 锯齿波 */
QWAVE_TYPE_ANTSAWTOOTH   /* 反向锯齿波 */
QWAVE_TYPE_NOISE         /* 伪随机噪声 */
QWAVE_TYPE_SQUARE        /* 方波 */
```

### 函数

#### `qwave_init`

```c
int qwave_init(QWaveCore *core, QWaveType type, fp_t fs,
               fp_t frq, fp_t bias, uint32_t seed);
```

初始化发生器。`fs` 是采样率（Hz），`frq` 是波形频率（Hz），`bias` 是直流偏置，`seed` 是噪声种子。成功返回 `0`；发生器指针为空、波形类型无效，或 `fs`/`frq` 不大于零时返回 `-1`。

初始化将振幅设为 `1`、相位和最近输出设为 `0`。若 `seed == 0`，内部使用固定的非零种子 `2463534242`，以避免 xorshift32 的零状态。

#### `qwave_output`

```c
fp_t qwave_output(QWaveCore *core);
```

生成并返回下一个输出样本，并更新 `gen->output` 和内部相位/随机状态。传入空指针时返回 `0`。该函数没有单独的错误码；正常输出也可能为零。

#### `qwave_signal_set`

```c
int qwave_signal_set(QWaveCore *core, QWaveType type);
```

切换波形类型。成功返回 `0`；空指针或无效波形类型返回 `-1`。切换时将相位和最近输出清零；不会重新初始化噪声 PRNG 状态。

#### `qwave_fs_set`

```c
int qwave_fs_set(QWaveCore *core, fp_t fs);
```

设置采样率并更新采样间隔 `ts`。成功返回 `0`；空指针或 `fs <= 0` 返回 `-1`。不会重置当前相位。

#### `qwave_frq_set`

```c
int qwave_frq_set(QWaveCore *core, fp_t frq);
```

设置频率并更新周期和半周期；当前相位会按新周期回绕。成功返回 `0`；空指针或 `frq <= 0` 返回 `-1`。

#### `qwave_amp_set`

```c
int qwave_amp_set(QWaveCore *core, fp_t amp);
```

设置波形振幅。成功返回 `0`；空指针或 `amp < 0` 返回 `-1`。零振幅有效，此时输出只剩偏置。

#### `qwave_bias_set`

```c
int qwave_bias_set(QWaveCore *core, fp_t bias);
```

设置直流偏置。成功返回 `0`；空指针返回 `-1`。

#### `qwave_reset`

```c
int qwave_reset(QWaveCore *core);
```

将相位和最近输出清零，并把 PRNG 状态恢复到初始化时对应的种子。成功返回 `0`；空指针返回 `-1`。它不会更改波形类型、采样率、频率、振幅或偏置。

## 按块生成样本

发生器可以在音频/控制循环中反复调用，也可以填充缓冲区：

```c
QWaveCore core;
fp_t samples[256];

if(qwave_init(&core, QWAVE_TYPE_TRIANGLE, 48000.0f, 440.0f, 0.0f, 1234) != 0) {
    return 1;
}

for(unsigned int i = 0; i < sizeof(samples) / sizeof(samples[0]); ++i) {
    samples[i] = qwave_output(&core);
}
```

## 注意事项

- `qwave_init()` 与 `qwave_fs_set()` 不强制检查 `fs > 5 * frq`；应由调用方选择合适采样率。
- 结构体字段是公开的。直接修改字段可能造成 `period`、`half_period` 或 `ts` 与配置不一致；频率、采样率等参数应使用对应 setter。
- `QWaveCore` 实例各自维护状态；同一个实例不要在多个线程中无同步地并发调用或修改。
- 伪随机噪声是确定性的，不是密码学随机数。相同种子和调用序列会产生相同噪声序列；`qwave_reset()` 会重新开始该序列。
- `qwave_output(NULL)` 返回 `0` 作为保护性行为，不能用返回值判断调用成功与否。

## 许可

本项目采用 MIT License，详见 [LICENSE](LICENSE)。
