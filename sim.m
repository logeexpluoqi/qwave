clc; close all;

frq = 1;
fs = 25;

T = 1/frq;
half_T = T/2;

ts = 1/fs;
t = 0 : ts : 2;

%% noise
close all;
ti = 0;
seed = int32(2463534242);
out = zeros(size(t));
for i = 1 : length(t) 
    x = seed;
    x = bitxor(bitshift(x, 13), x);
    x = bitxor(bitshift(x, -17), x);
    x = bitxor(bitshift(x, 5), x);
    seed = x;
    out(i) = (single(x) / 4294967295) * 2;
    ti = ti + ts;
    if ti >= T
        ti = ti - T;
    end
end

% 绘制信号
plot(t, out);
grid minor;
xlabel('时间');
ylabel('幅度');
title('信号仿真');

%% antsaw
close all;
ti = 0;
out = zeros(size(t));
for i = 1 : length(t) 
    norm = ti / T;
    out(i) = -norm;

    ti = ti + ts;
    if ti >= T
        ti = ti - T;
    end
end

% 绘制信号
plot(t, out);
grid minor;
xlabel('时间');
ylabel('幅度');
title('信号仿真');

%% saw
close all;
ti = 0;
out = zeros(size(t));
for i = 1 : length(t) 
    norm = ti / T;
    out(i) = norm;

    ti = ti + ts;
    if ti >= T
        ti = ti - T;
    end
end

% 绘制信号
plot(t, out);
grid minor;
xlabel('时间');
ylabel('幅度');
title('信号仿真');


%% tri
close all;
ti = 0;
out = zeros(size(t));
for i = 1 : length(t) 
    x = ti/T;
    if x < 0.25
        out(i) = 4 * x;
    elseif x < 0.75
        out(i) = 2 - 4 * x;
    else
        out(i) = -4 + 4 * x;
    end

    ti = ti + ts;
    if ti >= T
        ti = ti - T;
    end
end

% 绘制信号
plot(t, out);
grid minor;
xlabel('时间');
ylabel('幅度');
title('信号仿真');

%% sqr
close all;
ti = 0;
out = zeros(size(t));
for i = 1 : length(t) 
    
    if ti < half_T
        out(i) = -1;
    else
        out(i) = 1;
    end

    ti = ti + ts;
    if ti >= T
        ti = ti - T;
    end
end

% 绘制信号
plot(t, out);
grid minor;
xlabel('时间');
ylabel('幅度');
title('信号仿真');

%% sin
close all;
ti = 0;
out = zeros(size(t));
for i = 1 : length(t) 
    
    x = (ti / T) * 360;
    out(i) = fsin(x);

    ti = ti + ts;
    if ti >= T
        ti = ti - T;
    end
end

% 绘制信号
plot(t, out);
grid minor;
xlabel('时间');
ylabel('幅度');
title('信号仿真');
%% fsin
function out = fsin(t)
    fast_sin_table = [
    0.0, 0.017452406, 0.034899497, 0.052335956, 0.069756474, ...
    0.087155743, 0.104528463, 0.121869343, 0.139173101, 0.156434465,....
    0.173648178, 0.190809, 0.207911, 0.224951, 0.241922, ...
    0.258819, 0.275637, 0.292372, 0.309017, 0.325568, ...
    0.342020, 0.358368, 0.374607, 0.390731, 0.406737, ...
    0.422618, 0.438371, 0.453990, 0.469472, 0.48481, ...
    0.5, 0.515038, 0.529919, 0.544639, 0.559193, ...
    0.573576, 0.587785, 0.601815, 0.615661, 0.62932, ...
    0.642788, 0.656059, 0.669131, 0.681998, 0.694658, ...
    0.707107, 0.71934, 0.731354, 0.743145, 0.75471, ...
    0.766044, 0.777146, 0.788011, 0.798636, 0.809017, ...
    0.819152, 0.829038, 0.838671, 0.848048, 0.857167, ...
    0.866025, 0.87462, 0.882948, 0.891007, 0.898794, ...
    0.906308, 0.913545, 0.920505, 0.927184, 0.93358, ...
    0.939693, 0.945519, 0.951057, 0.956305, 0.961262, ...
    0.965926, 0.970296, 0.97437, 0.978148, 0.981627, ...
    0.984808, 0.987688, 0.990268, 0.992546, 0.994522, ...
    0.996195, 0.997564, 0.99863, 0.999391, 0.999848, ...
    1.0];


    t = fmodf(t, 360);
    
    sign = 1;
    % 使用sin的周期对称性
    if t > 180
        t = t - 180;
        sign = -1;
    end
    if t > 90
        t = 180 - t;
    end
    
    % 计算查找表索引和插值因子
    index = floor(t);  % 1°分辨率
    if index >= 90
        out = sign * fast_sin_table(91);
    else
        fraction = t - index;
        out = fast_sin_table(index + 1) * (1 - fraction) + fast_sin_table(index + 2) * fraction;
        out = sign * out;
    end
end

function result = fmodf(x, y)
    if y == 0
        result = NaN;
    else
        result = mod(x, y);
    end
end
