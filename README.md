InSAR & Groundwater Seasonal Deformation Analysis Toolpack

基于 InSAR 与地下水位的季节性形变模型与水文响应分析工具包

📖 English Version

Overview

This repository provides a modular, production-grade Python workflow for decomposing Interferometric Synthetic Aperture Radar (InSAR) time-series displacement and quantifying its hydrogeological coupling with groundwater-level dynamics.

The methodology is inspired by and aligned with research from the InSAR Geohazard Team at the China University of Mining and Technology (Beijing) (Tang et al., 2022, 2023, 2026). It enables automated extraction of linear subsidence rates, annual/semi-annual seasonal amplitudes and phases, cross-wavelet phase lags (Time Lag), and calibrated elastic skeletal storage coefficients ($S_{ke\_cal}$).

🌟 Key Features

6-Parameter Least-Squares Fitting (src/fitting.py):
Decomposes displacement time series into linear trends and periodic (annual + semi-annual) components. Includes rank-deficiency detection, RMSE computation, and phase unit conversions (radians/months).

Cross-Wavelet Transform (XWT) Analysis (src/wavelet_analysis.py):
Calculates cross-wavelet power spectra between surface deformation and groundwater head changes using Complex Morlet wavelets to extract phase lag ($\Delta t$) at dominant annual cycles.

Fault-Isolated Batch Processing (scripts/batch_process.py):
Processes multi-well CSV/DataFrame tables (long & wide formats) with automatic column mapping, time-lag sequence alignment, and per-well error isolation.

Data-Driven Hydrogeology:
Calibrates elastic skeletal storage coefficients ($S_{ke\_cal} = \frac{dH(t+\Delta t)}{dh}$) and estimates seasonal and long-term groundwater storage variations ($\Delta V$).

📐 Mathematical Formulation

1. Seasonal Deformation Decomposition Model

$$y(t) = y_0 + vt + \sum_{k=1}^{2} \left[ a_k \sin(2\pi k t) + b_k \cos(2\pi k t) \right] + \varepsilon(t)$$

Linear Rate ($v$): $\text{mm/year}$

Seasonal Amplitude ($A_k$): $A_k = \sqrt{a_k^2 + b_k^2}$

Seasonal Phase ($\phi_k$): $\phi_k = \arctan\left(\frac{a_k}{b_k}\right)$

2. Time-Lag Calibrated Storativity

$$S_{ke\_cal} = \frac{dH(t + \Delta t)}{dh}$$

Where $\Delta t$ is the poroelastic delay time derived from the annual phase angle $\Delta\Phi$ of the Cross-Wavelet Transform:

$$\Delta t = \frac{\Delta\Phi \cdot T}{2\pi}$$

📁 Project Structure

.
├── src/                        # Core algorithmic package
│   ├── __init__.py             # Public API exports
│   ├── fitting.py              # Least-squares seasonal decomposition
│   ├── wavelet_analysis.py     # CWT & XWT time-lag estimation
│   └── aquifer_params.py       # Ske_cal calibration and storage change
├── scripts/                    # Execution & batch-processing pipelines
│   ├── batch_process.py        # Multi-well batch runner with fault isolation
│   └── plot_results.py         # Spatial & scatter visualization scripts
├── data/                       # Input CSVs and synthetic datasets
├── results/                    # Output CSV summaries and exported figures
└── README.md                   # Project documentation



🛠️ Quick Start

1. Installation

Ensure Python 3.8+ is installed. Install dependencies in a virtual environment:

python3 -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate
pip install numpy scipy pywavelets pandas matplotlib geopandas



2. Run Single-Point / Module Test

import numpy as np
from src.fitting import fit_seasonal_model
from src.wavelet_analysis import compute_xwt_lag

# Synthetic test
t = np.arange(2017.0, 2024.0, 0.03)
y_insar = 5.0 - 30.0 * (t - 2017.0) + 15.0 * np.sin(2 * np.pi * t) + 10.0 * np.cos(2 * np.pi * t)
y_gw = 100.0 - 1.5 * (t - 2017.0) + 12.0 * np.cos(2 * np.pi * t - 1.0)

res = fit_seasonal_model(t, y_insar)
print(f"Velocity: {res.v:.2f} mm/yr, Amplitude: {res.amp_annual:.2f} mm")

xwt_res = compute_xwt_lag(y_insar, y_gw, dt=0.03)
print(f"Time Lag Δt: {xwt_res.time_lag_days:.1f} days")



3. Execute Batch Pipeline

python scripts/batch_process.py --input data/wells_timeseries.csv --output results/multi_wells_summary.csv



🇨🇳 中文版 (Chinese Version)

项目简介

本仓库提供了一套模块化、生产级的 Python 数据分析流水线，旨在剥离 InSAR 形变时间序列中的长期趋势与季节性信号，并深入量化地表形变与地下水位动态之间的水文地质耦合关系。

本项目依据中国矿业大学（北京）InSAR 地质灾害研究团队（Tang et al., 2022, 2023, 2026）提出的学术框架构建，实现了线性沉降速率、年/半年周期振幅与相位、交叉小波响应延迟天数（Time Lag）以及校正后弹性骨架贮水系数（$S_{ke\_cal}$）的自动化提取。

🌟 核心特性

六参数最小二乘拟合 (src/fitting.py)：
将形变时序拆解为线性趋势与周期分量（年+半年），内置秩亏检测、RMSE 计算及月度相位单位转换。

交叉小波分析 XWT (src/wavelet_analysis.py)：
基于复莫莱（Complex Morlet）小波计算形变与水位的交叉功率谱，精准计算 1 年主周期上的响应延迟天数 ($\Delta t$)。

故障隔离批处理流水线 (scripts/batch_process.py)：
支持多监测井 CSV/DataFrame 长宽表输入，具备自动列名识别、时序延时对齐及单井错误隔离能力。

水文参数定量评估：
基于延时校正公式 $S_{ke\_cal} = \frac{dH(t+\Delta t)}{dh}$ 评估骨架贮水系数，进而估算区域地下水季节性与长期储量变化量 ($\Delta V$)。

📐 数学原理

1. 季节性形变分解模型

$$y(t) = y_0 + vt + \sum_{k=1}^{2} \left[ a_k \sin(2\pi k t) + b_k \cos(2\pi k t) \right] + \varepsilon(t)$$

线性速率 ($v$): 单位 $\text{mm/年}$

周期振幅 ($A_k$): $A_k = \sqrt{a_k^2 + b_k^2}$

周期相位 ($\phi_k$): $\phi_k = \arctan\left(\frac{a_k}{b_k}\right)$

2. 滞后校正贮水系数

$$S_{ke\_cal} = \frac{dH(t + \Delta t)}{dh}$$

其中响应延迟天数 $\Delta t$ 来源于交叉小波分析得到的年周期相角差 $\Delta\Phi$：

$$\Delta t = \frac{\Delta\Phi \cdot T}{2\pi}$$

📁 项目结构

.
├── src/                        # 核心算法模块包
│   ├── __init__.py             # 公共接口导出
│   ├── fitting.py              # 最小二乘季节拟合模型
│   ├── wavelet_analysis.py     # 连续与交叉小波(CWT/XWT)延迟计算
│   └── aquifer_params.py       # 延时校正 Ske 与储量估算
├── scripts/                    # 批处理与可视化脚本
│   ├── batch_process.py        # 支持单井隔离的多点批处理管道
│   └── plot_results.py         # 空间热力图与散点关系绘制
├── data/                       # 输入 CSV 格式及合成测试数据
├── results/                    # 输出的汇总 CSV 表格与导出图表
└── README.md                   # 中英文项目文档



🛠️ 快速上手

1. 环境配置

使用 Python 3.8+ 环境，建议在虚拟环境中安装依赖：

python3 -m venv venv
source venv/bin/activate  # Windows 运行: venv\Scripts\activate
pip install numpy scipy pywavelets pandas matplotlib geopandas



2. 批量处理流水线运行

python scripts/batch_process.py --input data/wells_timeseries.csv --output results/multi_wells_summary.csv



📚 References & Acknowledgments / 参考文献与致谢

Tang, W., Zhao, X., Liu, Y., Liu, R., Fan, B., Zhu, T. (2026). Seasonal and Long-Term Aquifer Deformation Revealed by InSAR Remote Sensing: Implications for Groundwater Recharge in an Arid Oasis of the Northern Tianshan Piedmont.

Tang, W., Zhao, X., Bi, G., et al. (2023). Quantifying seasonal ground deformation in Taiyuan basin, China, by Sentinel-1 InSAR time series analysis. Journal of Hydrology, 622, 129654.

Tang, W., Zhao, X., Motagh, M., et al. (2022). Land subsidence and rebound in the Taiyuan basin, northern China, in the context of inter-basin water transfer and groundwater management. Remote Sensing of Environment, 269, 112792.

InSAR Geohazard Team, College of Geoscience and Surveying Engineering, China University of Mining and Technology (Beijing).

Maintained for the Undergraduate Innovation and Entrepreneurship Training Program (大创项目).
